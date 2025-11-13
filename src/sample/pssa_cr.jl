
# pssa_cr.jl
#
# Exact continuous-time sampler using composition–rejection over groups,
# in the spirit of Partial-Propensity SSA with Composition–Rejection (PSSA-CR).
#
# Design goals:
# - Match the sampler interface shown in firsttofire.jl (enable!, next, fire!, disable!, reset!, clone, copy_clocks!, jitter!, haskey/keys/length/getindex).
# - Work with Exponential clocks (Markov jump processes). Reject non-Exponential.
# - Provide exact sampling of the next reaction and time: Δt ~ Exp(a0), reaction chosen via composition–rejection within groups.
# - Be usable without reaction-network structure; optional manual grouping API is provided to recover PSSA-style efficiency when you have owner/group metadata.
#
# SIMPLIFICATIONS FROM THE ORIGINAL PSSA-CR ALGORITHM (Ramaswamy & Sbalzarini 2010):
#
# This implementation is a SIMPLIFIED, GENERAL-PURPOSE variant of PSSA-CR that trades some theoretical
# performance guarantees for broader applicability:
#
# 1. NO TRUE PARTIAL PROPENSITIES:
#    - Paper: Groups reactions by reactant species, factoring propensities as πμ(i) = nᵢ × cμ
#      (partial propensity = population × rate constant)
#    - Here: Groups generic "clocks" by hash value, stores full rates λ (not factored)
#    - Impact: Works for any exponential process, not just chemical reactions with species counts
#
# 2. NO DYADIC BINNING:
#    - Paper: Uses log₂(λmax/λmin)+1 bins per group for O(1) composition-rejection
#    - Here: Uses LINEAR SEARCH over groups (lines 257-270)
#    - Impact: Composition step is O(G) instead of O(log G), where G = number of groups
#    - Trade-off: Simpler implementation, still efficient for moderate G (default 64)
#
# 3. NO DEPENDENCY GRAPH OVER SPECIES:
#    - Paper: Maintains species dependency graph to update O(1) partial propensities for weakly-coupled networks
#    - Here: No dependency tracking; each enable!/disable! only updates affected group
#    - Impact: Update cost is O(1) per operation (not O(N) for strongly-coupled networks)
#
# 4. HASH-BASED GROUPING (NOT REACTION STRUCTURE):
#    - Paper: Groups reactions by their reactant species (exploiting chemical network structure)
#    - Here: Groups clocks by hash(clock_id) % ngroups
#    - Impact: Grouping is general-purpose but may have worse rejection rates than structure-aware grouping
#    - Mitigation: User can manually assign groups via assign_group!() if structure is known
#
# 5. THEORETICAL COMPLEXITY:
#    - Paper PSSA-CR: O(1) per step for weakly-coupled, O(N) for strongly-coupled networks
#    - This implementation: O(G) per step where G = ngroups (default 64)
#    - Note: Since G is typically constant and small, this is effectively O(1) in practice
#
# WHEN TO USE THIS IMPLEMENTATION:
# - ✓ General exponential processes (not just chemical reactions)
# - ✓ Moderate number of groups (G ≤ 100 works well)
# - ✓ When simplicity and code clarity matter
# - ✗ Extremely large-scale chemical networks where true O(1) scaling is critical
# - ✗ When you need absolute minimum rejection rate (use structure-aware grouping)
#
# Notes:
# - This sampler keeps no per-clock scheduled times. It caches *only one* "next event" (time,clock)
#   to make `next` idempotent until `fire!`, `enable!`, `disable!`, or `jitter!` invalidates it.
# - Grouping: by default, clocks are hashed into a fixed number of groups. You can override the group
#   assignment before enabling a clock via `assign_group!` for better acceptance behavior.
#
# References:
# R. Ramaswamy & I. F. Sbalzarini, "A partial-propensity variant of the composition-rejection SSA", J. Chem. Phys. 132, 044102 (2010).
# A. Slepoy, A. P. Thompson, S. J. Plimpton, "A constant-time kinetic Monte Carlo algorithm...", J. Chem. Phys. 128, 205101 (2008).
#

export PSSACR

"""
    PSSACR{KeyType,TimeType}(; ngroups::Int=64)

Exact, continuous-time sampler using composition–rejection over groups.
Meant to be a drop-in replacement where `FirstToFire` is used.

Assumptions:
- Each enabled clock is Exponential with rate λ (i.e., distribution `Exponential(θ)` where λ=1/θ).
- The sampler maintains group-wise sums and maxima of rates.
- Time to next event is drawn from `Exp(∑λ)`, and the firing key is drawn by a two-stage
  selection: group by probability ∝ group-sum, then within-group by rejection from a uniform
  proposal with acceptance `λ/λ_max_group`.

Interface compatibility:
- `next` is idempotent: it caches one upcoming `(time,clock)` until invalidated by `fire!`,
  `enable!`, `disable!`, or `jitter!`.

Performance notes:
- Choose `ngroups` so that groups remain small and rate magnitudes similar.
  If you possess "owner" metadata (partial-propensity style), call `assign_group!`
  *before* `enable!` to place clocks with common owners into the same group.
"""
mutable struct PSSACR{K,T} <: SSA{K,T}
    # Per-clock current rate λ_k
    rates::Dict{K,T}
    # Group membership and positions
    group_of::Dict{K,Int}            # clock -> group id in 1:ngroups
    pos_in_group::Dict{K,Int}        # clock -> index inside groups[group]
    groups::Vector{Vector{K}}        # groups[g] is a Vector of keys in group g
    # Group statistics
    group_sum::Vector{T}             # Σ_k∈g λ_k
    group_max::Vector{T}             # max_k∈g λ_k (for rejection)
    total_rate::T                    # Σ_g group_sum[g]
    # Cached next event for idempotent `next`
    cached_next::Union{Nothing,OrderedSample{K,T}}
end

# Constructor
function PSSACR{K,T}(; ngroups::Int=64) where {K,T}
    groups = [Vector{K}() for _ in 1:ngroups]
    zerosT = zero(T)
    PSSACR{K,T}(
        Dict{K,T}(),
        Dict{K,Int}(),
        Dict{K,Int}(),
        groups,
        fill(zerosT, ngroups),
        fill(zerosT, ngroups),
        zerosT,
        nothing,
    )
end

# Shallow clone of structure (no clocks), preserving ngroups
function clone(s::PSSACR{K,T}) where {K,T}
    PSSACR{K,T}(ngroups=length(s.groups))
end

# Reset all internal state
function reset!(s::PSSACR{K,T}) where {K,T}
    empty!(s.rates)
    empty!(s.group_of)
    empty!(s.pos_in_group)
    for g in eachindex(s.groups)
        empty!(s.groups[g])
        s.group_sum[g] = zero(T)
        s.group_max[g] = zero(T)
    end
    s.total_rate = zero(T)
    s.cached_next = nothing
    return s
end

# Deep copy clocks and statistics from src into dst
function copy_clocks!(dst::PSSACR{K,T}, src::PSSACR{K,T}) where {K,T}
    # Validate matching structure
    if length(dst.groups) != length(src.groups)
        throw(ArgumentError("copy_clocks!: destination has $(length(dst.groups)) groups but source has $(length(src.groups))"))
    end

    dst.rates = copy(src.rates)
    dst.group_of = copy(src.group_of)
    dst.pos_in_group = copy(src.pos_in_group)
    dst.groups = [copy(v) for v in src.groups]
    dst.group_sum = copy(src.group_sum)
    dst.group_max = copy(src.group_max)
    dst.total_rate = src.total_rate
    dst.cached_next = src.cached_next === nothing ? nothing :
        OrderedSample{K,T}(getfield(src.cached_next, :clock), getfield(src.cached_next, :time))
    return dst
end

# Optional helper: set a group for a clock BEFORE enabling it.
function assign_group!(s::PSSACR{K,T}, clock::K, group::Int) where {K,T}
    @boundscheck 1 <= group <= length(s.groups) || throw(ArgumentError("group out of range"))
    s.group_of[clock] = group
    return s
end

# Internal: invalidate cached next event
@inline function _invalidate!(s::PSSACR)
    s.cached_next = nothing
end

# Internal: recompute max rate in a group (O(group size))
function _recompute_group_max!(s::PSSACR{K,T}, g::Int) where {K,T}
    maxv = zero(T)
    for k in s.groups[g]
        rk = s.rates[k]
        if rk > maxv
            maxv = rk
        end
    end
    s.group_max[g] = maxv
end

# Internal: add a brand-new key with rate λ to its group
function _insert_new!(s::PSSACR{K,T}, clock::K, λ::T) where {K,T}
    g = get(s.group_of, clock, 1 + (abs(hash(clock)) % length(s.groups)))
    push!(s.groups[g], clock)
    s.pos_in_group[clock] = length(s.groups[g])
    s.group_of[clock] = g
    s.rates[clock] = λ
    s.group_sum[g] += λ
    if λ > s.group_max[g]; s.group_max[g] = λ; end
    s.total_rate += λ
end

# Internal: update existing key's rate
function _update_rate!(s::PSSACR{K,T}, clock::K, newλ::T) where {K,T}
    oldλ = s.rates[clock]
    if oldλ == newλ; return; end
    g = s.group_of[clock]
    s.rates[clock] = newλ
    s.group_sum[g] += (newλ - oldλ)
    s.total_rate += (newλ - oldλ)
    if newλ > s.group_max[g]
        s.group_max[g] = newλ
    elseif oldλ == s.group_max[g] && newλ < oldλ
        _recompute_group_max!(s, g)
    end
end

# Required interface: jitter! → just invalidate the cached sample.
function jitter!(s::PSSACR{K,T}, when::T, rng::AbstractRNG) where {K,T}
    _invalidate!(s)
end

# Required interface: enabling or updating a clock.
function enable!(s::PSSACR{K,T},
                 clock::K,
                 distribution::UnivariateDistribution,
                 te::T,
                 when::T,
                 rng::AbstractRNG) where {K,T}
    if !(distribution isa Exponential)
        throw(ArgumentError("PSSACR only supports Exponential distributions (got $(typeof(distribution)))."))
    end

    # Distributions.jl: Exponential(θ) has rate 1/θ
    λ = convert(T, rate(distribution))

    if haskey(s.rates, clock)
        _update_rate!(s, clock, λ)
    else
        _insert_new!(s, clock, λ)
    end
    _invalidate!(s)
    return s
end

# Required interface: disable a clock completely.
function disable!(s::PSSACR{K,T}, clock::K, when::T) where {K,T}
    haskey(s.rates, clock) || throw(KeyError(clock))
    g = s.group_of[clock]
    pos = s.pos_in_group[clock]
    oldλ = s.rates[clock]

    # swap-remove within group
    lastidx = length(s.groups[g])
    if pos != lastidx
        klast = s.groups[g][lastidx]
        s.groups[g][pos] = klast
        s.pos_in_group[klast] = pos
    end
    pop!(s.groups[g])
    delete!(s.pos_in_group, clock)
    delete!(s.group_of, clock)

    # update stats
    s.group_sum[g] -= oldλ
    s.total_rate -= oldλ
    if oldλ == s.group_max[g]
        _recompute_group_max!(s, g)
    end
    delete!(s.rates, clock)

    _invalidate!(s)
    return s
end

# Required interface: `fire!` has no internal scheduled time to remove.
# We simply invalidate the cached next event.
function fire!(s::PSSACR{K,T}, clock::K, when::T) where {K,T}
    _invalidate!(s)
    return s
end

# Required interface: compute (time, clock) of the next event WITHOUT removing it.
function next(s::PSSACR{K,T}, when::T, rng::AbstractRNG) where {K,T}
    # If cached, return it.
    if s.cached_next !== nothing
        t = getfield(s.cached_next, :time)
        k = getfield(s.cached_next, :clock)
        return (t, k)
    end

    # No active clocks
    if isempty(s.rates) || s.total_rate <= zero(T)
        return (typemax(T), nothing)
    end

    # Time increment
    Δt = rand(rng, Exponential(inv(s.total_rate)))
    tfire = when + Δt

    # Select group ∝ group_sum
    target = rand(rng) * s.total_rate
    acc = zero(T)
    gsel = 0
    @inbounds for g in 1:length(s.groups)
        gs = s.group_sum[g]
        if gs <= zero(T); continue; end
        acc += gs
        if target <= acc
            gsel = g
            break
        end
    end

    # Safety: if roundoff picked an empty or zero-sum group, find a non-empty one.
    if gsel == 0
        for g in 1:length(s.groups)
            if s.group_sum[g] > zero(T) && !isempty(s.groups[g])
                gsel = g
                break
            end
        end
        if gsel == 0
            return (typemax(T), nothing) # nothing enabled effectively
        end
    end

    # Within-group selection via rejection from uniform proposal.
    # Proposal: pick a key uniformly from groups[gsel]. Accept with prob λ/λ_max.
    λmax = s.group_max[gsel]
    @assert λmax >= zero(T)
    chosen_key = nothing
    keys_in_group = s.groups[gsel]
    @inbounds while true
        idx = rand(rng, 1:length(keys_in_group))
        k = keys_in_group[idx]
        λ = s.rates[k]
        if λ > zero(T) && (rand(rng) * λmax) <= λ
            chosen_key = k
            break
        end
        # else retry
    end

    ev = OrderedSample{K,T}(chosen_key, tfire)
    s.cached_next = ev
    return (tfire, chosen_key)
end

# For compatibility with the example sampler: expose stored time for a given key.
# Here we only cache the *single* next event. We return its time iff it belongs to `clock`.
function Base.getindex(s::PSSACR{K,T}, clock::K) where {K,T}
    s.cached_next === nothing && throw(KeyError(clock))
    key = getfield(s.cached_next, :clock)
    if key == clock
        return getfield(s.cached_next, :time)
    else
        throw(KeyError(clock))
    end
end

# Housekeeping
function Base.keys(s::PSSACR{K,T}) where {K,T}
    return collect(Base.keys(s.rates))
end

Base.length(s::PSSACR{K,T}) where {K,T} = length(s.rates)
Base.haskey(s::PSSACR{K,T}, clock::K) where {K,T} = haskey(s.rates, clock)
Base.haskey(::PSSACR{K,T}, ::Any) where {K,T} = false

enabled(s::PSSACR) = collect(Base.keys(s.rates))
isenabled(s::PSSACR{K,T}, clock::K) where {K,T} = haskey(s.rates, clock)
