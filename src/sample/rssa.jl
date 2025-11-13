
# rssa.jl
#
# Exact Rejection-based Stochastic Simulation Algorithm (RSSA)
# for continuous-time Markov jump processes (exponential clocks).
#
# It matches the sampler interface shown in firsttofire.jl:
#   next, enable!, fire!, disable!, reset!, clone, copy_clocks!,
#   jitter!, haskey/keys/length/getindex.
#
# Algorithmic core:
# - Maintain per-clock *true* rate a_i and a certified upper bound \bar a_i >= a_i.
# - Maintain Ā = sum_i \bar a_i and a Fenwick tree over { \bar a_i } for O(log N) sampling.
# - Draw candidate times from Exp(Ā). Select candidate clock i by Categorical(\bar a_i/Ā).
# - Accept with probability a_i / \bar a_i; otherwise reject and continue thinning.
#
# Exactness: standard thinning of a Poisson process with rate Ā, with acceptance a_i/\bar a_i,
# yields the target Markov jump process (homogeneous propensities). See Thanh et al. (2014, 2015).
#
# Notes:
# - This implementation targets time-homogeneous propensities (Exponential only).
#   For time-dependent rates (tRSSA), you need piecewise-time envelopes and an integral sampler;
#   these can be added as a thin extension without changing the public interface.
#
using Random
using Distributions: UnivariateDistribution, Exponential, rate

export RSSA, set_bound!, set_global_bound_factor!

"""
    RSSA{KeyType,TimeType}(; bound_factor=1.05)

Rejection-based SSA with global Fenwick tree for candidate selection.
- `bound_factor` ≥ 1.0 controls default upper bounds: \\bar a_i ← max(\\bar a_i, bound_factor * a_i).
  Set to 1.0 for no rejections (reduces to direct-method timing with tree selection).
"""
mutable struct RSSA{K,T} <: SSA{K,T}
    idx_of::Dict{K,Int}           # key → index (stable; indices are never reused)
    keys_vec::Vector{K}           # index → key
    present::BitVector            # index enabled?
    a::Vector{T}                  # true rates a_i
    abar::Vector{T}               # bounds \bar a_i ≥ a_i
    bit::Vector{T}                # Fenwick tree over abar
    Abar::T                       # sum(abar) of enabled indices
    cached_next::Union{Nothing,OrderedSample{K,T}}
    bound_factor::T               # default multiplicative bound factor (≥1)
end

function RSSA{K,T}(; bound_factor=1.05) where {K,T}
    bf = convert(T, bound_factor)
    bf < one(T) && (bf = one(T))
    RSSA{K,T}(
        Dict{K,Int}(),
        K[],
        BitVector(),
        T[],
        T[],
        T[],
        zero(T),
        nothing,
        bf,
    )
end

# Clone with empty state
clone(s::RSSA{K,T}) where {K,T} = RSSA{K,T}(; bound_factor=s.bound_factor)

# Reset everything
function reset!(s::RSSA{K,T}) where {K,T}
    empty!(s.idx_of)
    empty!(s.keys_vec)
    empty!(s.present)
    empty!(s.a)
    empty!(s.abar)
    empty!(s.bit)
    s.Abar = zero(T)
    s.cached_next = nothing
    return s
end

# Deep copy clocks and internal structures
function copy_clocks!(dst::RSSA{K,T}, src::RSSA{K,T}) where {K,T}
    dst.idx_of = copy(src.idx_of)
    dst.keys_vec = copy(src.keys_vec)
    dst.present = copy(src.present)
    dst.a = copy(src.a)
    dst.abar = copy(src.abar)
    dst.bit = copy(src.bit)
    dst.Abar = src.Abar
    dst.cached_next = src.cached_next === nothing ? nothing :
        OrderedSample{K,T}(getfield(src.cached_next, :key), getfield(src.cached_next, :time))
    dst.bound_factor = src.bound_factor
    return dst
end

# ---- Fenwick tree utilities over abar ----
@inline function _bit_add!(bit::Vector{T}, i::Int, δ::T) where {T}
    n = length(bit)
    while i <= n
        @inbounds bit[i] += δ
        i += i & -i
    end
end

@inline function _bit_sum(bit::Vector{T}, i::Int) where {T}
    s = zero(eltype(bit))
    while i > 0
        @inbounds s += bit[i]
        i -= i & -i
    end
    return s
end

# Find smallest index j with prefix sum ≥ target, assuming target ∈ (0, sum]
# Returns a value in [1, n]. If target exceeds the tree sum (due to floating-point drift),
# returns n to avoid out-of-bounds access.
function _bit_find(bit::Vector{T}, target::T) where {T}
    n = length(bit)
    idx = 0
    # largest power of two ≤ n
    bitmask = one(Int)
    while (bitmask << 1) <= n
        bitmask <<= 1
    end
    while bitmask != 0
        nxt = idx + bitmask
        if nxt <= n && bit[nxt] < target
            target -= bit[nxt]
            idx = nxt
        end
        bitmask >>= 1
    end
    # Clamp to valid range [1, n] to handle floating-point drift in Abar
    return min(idx + 1, n)
end

# ---- internal helpers ----
@inline function _invalidate!(s::RSSA)
    s.cached_next = nothing
end

# Ensure an index exists for `key`. Indices are append-only; removal just disables.
function _ensure_index!(s::RSSA{K,T}, key::K) where {K,T}
    idx = get(s.idx_of, key, 0)
    if idx != 0
        return idx
    end
    push!(s.keys_vec, key)
    push!(s.present, false)
    push!(s.a, zero(T))
    push!(s.abar, zero(T))
    push!(s.bit, zero(T))
    idx = length(s.keys_vec)
    s.idx_of[key] = idx
    return idx
end

# Public: set a tighter (or looser) *bound* for a clock. Ensures \bar a ≥ a.
function set_bound!(s::RSSA{K,T}, key::K, abar_new::T) where {K,T}
    idx = get(s.idx_of, key, 0)
    idx == 0 && throw(KeyError(key))
    # Enforce invariant: bound must be at least as large as true rate
    abar_new = max(abar_new, s.a[idx])
    abar_new < zero(T) && (abar_new = zero(T))
    δ = abar_new - s.abar[idx]
    s.abar[idx] = abar_new
    if s.present[idx]
        _bit_add!(s.bit, idx, δ)
        s.Abar += δ
    end
    _invalidate!(s)
    return s
end

# Public: set a new global bound factor and recompute all bounds for enabled clocks.
function set_global_bound_factor!(s::RSSA{K,T}, bf) where {K,T}
    s.bound_factor = convert(T, bf)
    s.bound_factor < one(T) && (s.bound_factor = one(T))
    # rebuild BIT and Abar
    fill!(s.bit, zero(T))
    s.Abar = zero(T)
    for idx in 1:length(s.keys_vec)
        if s.present[idx]
            s.abar[idx] = max(s.bound_factor * s.a[idx], eps(T))
            _bit_add!(s.bit, idx, s.abar[idx])
            s.Abar += s.abar[idx]
        end
    end
    _invalidate!(s)
    return s
end

# ---- interface methods ----

# No scheduled times to perturb; just drop cached sample.
function jitter!(s::RSSA{K,T}, when::T, rng::AbstractRNG) where {K,T}
    _invalidate!(s)
end

# Enable or update a clock. Only Exponential is supported (homogeneous).
function enable!(s::RSSA{K,T},
                 key::K,
                 distribution::UnivariateDistribution,
                 te::T,
                 when::T,
                 rng::AbstractRNG) where {K,T}
    if distribution isa Exponential
        λ = convert(T, rate(distribution))  # Distributions.jl: rate = 1/θ
    else
        throw(ArgumentError("RSSA only supports Exponential distributions (got $(typeof(distribution)))."))
    end

    idx = _ensure_index!(s, key)
    old_enabled = s.present[idx]
    oldabar = s.abar[idx]

    # Update true rate
    s.a[idx] = λ

    # Choose a default bound if needed
    newabar = max(oldabar, s.bound_factor * λ)
    if !old_enabled
        # enable
        s.present[idx] = true
        # Don't artificially inflate zero rates - allow abar to be zero
        # This prevents infinite loops when all enabled clocks have zero rates
        s.abar[idx] = newabar
        s.Abar += newabar
        _bit_add!(s.bit, idx, newabar)
    else
        # already enabled: adjust bound minimally to keep abar ≥ a
        if newabar < s.a[idx]
            newabar = s.a[idx]
        end
        if newabar != oldabar
            δ = newabar - oldabar
            s.abar[idx] = newabar
            s.Abar += δ
            _bit_add!(s.bit, idx, δ)
        end
    end
    _invalidate!(s)
    return s
end

# Disable a clock
function disable!(s::RSSA{K,T}, key::K, when::T) where {K,T}
    idx = get(s.idx_of, key, 0)
    idx == 0 && throw(KeyError(key))
    s.present[idx] || throw(KeyError(key))
    # remove its bound mass
    w = s.abar[idx]
    if w > zero(T)
        _bit_add!(s.bit, idx, -w)
        s.Abar -= w
    end
    s.present[idx] = false
    s.a[idx] = zero(T)
    s.abar[idx] = zero(T)
    _invalidate!(s)
    return s
end

# After firing, nothing to remove; just invalidate cached sample.
function fire!(s::RSSA{K,T}, key::K, when::T) where {K,T}
    _invalidate!(s)
    return s
end

# Idempotent: caches a single next (time,key) until invalidated by fire!/enable!/disable!/jitter!
function next(s::RSSA{K,T}, when::T, rng::AbstractRNG) where {K,T}
    if s.cached_next !== nothing
        ev = s.cached_next
        return (getfield(ev, :time), getfield(ev, :key))
    end

    if s.Abar <= zero(T)
        return (typemax(T), nothing)
    end

    t = when
    while true
        # candidate time from Exp(Abar)
        Δ = rand(rng, Exponential(inv(s.Abar)))
        t += Δ

        # pick candidate index by weights \bar a_i using Fenwick search
        u = rand(rng) * s.Abar
        j = _bit_find(s.bit, u)

        # in case of numerical corner cases, resample
        if j < 1 || j > length(s.keys_vec) || !s.present[j] || s.abar[j] <= zero(T)
            continue
        end

        # accept with probability a_j / \bar a_j
        aj = s.a[j]
        if aj <= zero(T)
            # zero true rate → always reject; continue thinning
            continue
        end

        if rand(rng) <= aj / s.abar[j]
            key = s.keys_vec[j]
            ev = OrderedSample{K,T}(key, t)
            s.cached_next = ev
            return (t, key)
        end
        # else: rejected virtual reaction ⇒ continue loop
    end
end

# For this sampler, we only cache the *single* next event.
# Return its time only if it belongs to `key`.
function Base.getindex(s::RSSA{K,T}, key::K) where {K,T}
    s.cached_next === nothing && throw(KeyError(key))
    k = getfield(s.cached_next, :key)
    k == key || throw(KeyError(key))
    return getfield(s.cached_next, :time)
end

# Introspection helpers
function Base.keys(s::RSSA{K,T}) where {K,T}
    out = Vector{K}()
    for (i,k) in enumerate(s.keys_vec)
        if s.present[i]
            push!(out, k)
        end
    end
    return out
end

Base.length(s::RSSA{K,T}) where {K,T} = count(s.present)
Base.haskey(s::RSSA{K,T}, key::K) where {K,T} = (get(s.idx_of, key, 0) |> i->(i!=0 && s.present[i]))
Base.haskey(::RSSA{K,T}, ::Any) where {K,T} = false

enabled(s::RSSA) = keys(s)
isenabled(s::RSSA{K,T}, key::K) where {K,T} = haskey(s, key)
