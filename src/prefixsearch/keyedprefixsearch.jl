import Base: in, setindex!, delete!
using Random


abstract type KeyedPrefixSearch end

"""
This decorator turns a Prefix Search algorithm into one that works for
arbitrary keys. This version only adds entries, so disabling a
clock sets its hazard to zero without removing it. If a simulation
re-enables the same set of clocks, this is the faster choice.
"""
struct KeyedKeepPrefixSearch{T,P} <: KeyedPrefixSearch
    # Map from clock name to index in propensity array.
    index::Dict{T, Int}
    # Map from index in propensity array to clock name.
    key::Vector{T}
    prefix::P
    function KeyedKeepPrefixSearch{T,P}(prefix::P) where {T,P}
        new{T,P}(Dict{T,Int}(), Vector{T}(), prefix)
    end
end


Base.length(kp::KeyedKeepPrefixSearch) = length(kp.index)

function Base.setindex!(kp::KeyedKeepPrefixSearch, val, clock)
    idx = get(kp.index, clock, 0)
    if idx != 0
        kp.prefix[idx] = val
    else
        push!(kp.prefix, val)
        kp.index[clock] = length(kp.prefix)
        push!(kp.key, clock)
        @assert length(kp.key) == length(kp.prefix)
    end
end


Base.delete!(kp::KeyedKeepPrefixSearch, clock) = kp.prefix[kp.index[clock]] = zero(Float64)
function Base.sum!(kp::KeyedKeepPrefixSearch)
    (length(kp.index) > 0) ? sum!(kp.prefix) : zero(Float64)
end


function choose(kp::KeyedKeepPrefixSearch, value)
    idx, hazard = choose(kp.prefix, value)
    return (kp.key[idx], hazard)
end


function Random.rand(
    rng::AbstractRNG, d::Random.SamplerTrivial{KeyedKeepPrefixSearch{T,P}}
    ) where {T,P}
    total = sum!(d[])
    choose(d[], rand(rng, Uniform{Float64}(0, total)))
end


"""
This decorator turns a Prefix Search algorithm into one that works for
arbitrary keys. This version reuses entries in the prefix search
after their clocks have been disabled. If the simulation moves through
a large key space, this will use less memory.
"""
struct KeyedRemovalPrefixSearch{T,P} <: KeyedPrefixSearch
    # Map from clock name to index in propensity array.
    index::Dict{T, Int}
    # Map from index in propensity array to clock name.
    key::Vector{T}
    free::Set{Int}
    prefix::P
    function KeyedRemovalPrefixSearch{T,P}(prefix::P) where {T,P}
        new{T,P}(Dict{T,Int64}(), Vector{T}(), Set{Int}(), prefix)
    end
end


Base.length(kp::KeyedRemovalPrefixSearch) = length(kp.index)

function Base.setindex!(kp::KeyedRemovalPrefixSearch, val, clock)
    idx = get(kp.index, clock, 0)
    if idx != 0
        kp.prefix[idx] = val
    elseif !isempty(kp.free)
        idx = pop!(kp.free)
        kp.index[clock] = idx
        kp.key[idx] = clock
        kp.prefix[idx] = val
    else
        push!(kp.prefix, val)
        kp.index[clock] = length(kp.prefix)
        push!(kp.key, clock)
        @assert length(kp.key) == length(kp.prefix)
    end
end


function Base.delete!(kp::KeyedRemovalPrefixSearch, clock)
    idx = kp.index[clock]
    kp.prefix[idx] = zero(Float64)
    delete!(kp.index, clock)
    # kp.key[idx] is now out of date.
    push!(kp.free, idx)
end


function Base.sum!(kp::KeyedRemovalPrefixSearch)
    (length(kp.index) > 0) ? sum!(kp.prefix) : zero(Float64)
end

function choose(kp::KeyedRemovalPrefixSearch, value)
    idx, hazard = choose(kp.prefix, value)
    return (kp.key[idx], hazard)
end


function Random.rand(
    rng::AbstractRNG, d::Random.SamplerTrivial{KeyedRemovalPrefixSearch{T,P}}
    ) where {T,P}
    total = sum!(d[])
    choose(d[], rand(rng, Uniform{Float64}(0, total)))
end
