import Base: in, setindex!, delete!
using Random


"""
This decorator turns a Prefix Search algorithm into one that works for
arbitrary keys. This version only adds entries, so disabling a
clock sets its hazard to zero without removing it. If a simulation
re-enables the same set of clocks, this is the faster choice.
"""
struct KeyedPrefixSearch{T,P,F<:Function}
    # Map from clock name to index in propensity array.
    index::Dict{T, Int}
    # Map from index in propensity array to clock name.
    key::Vector{T}
    key_space::F
    prefix::P
    function KeyedPrefixSearch{T,P}(prefix::P) where {T,P}
        yes = _ -> true
        new{T,P,typeof(yes)}(Dict{T,Int}(), Vector{T}(), yes, prefix)
    end
end


Base.in(kp::KeyedPrefixSearch, clock_id) = kp.key_space(clock_id)


function Base.setindex!(kp::KeyedPrefixSearch, val, clock)
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


Base.delete!(kp::KeyedPrefixSearch, clock) = kp.prefix[kp.index[clock]] = zero(Float64)
function Base.sum!(kp::KeyedPrefixSearch)
    (length(kp.index) > 0) ? sum!(kp.prefix) : zero(Float64)
end


function choose(kp::KeyedPrefixSearch, value)
    idx, hazard = choose(kp.prefix, value)
    return (kp.key[idx], hazard)
end


function Random.rand(
    rng::AbstractRNG, d::Random.SamplerTrivial{KeyedPrefixSearch{T,P,F}}
    ) where {T,P,F}
    total = sum!(d[])
    choose(d[], rand(rng, Uniform{Float64}(0, total)))
end


"""
This decorator turns a Prefix Search algorithm into one that works for
arbitrary keys. This version reuses entries in the prefix search
after their clocks have been disabled. If the simulation moves through
a large key space, this will use less memory.
"""
struct KeyedRemovalPrefixSearch{T,P,F}
    # Map from clock name to index in propensity array.
    index::Dict{T, Int}
    # Map from index in propensity array to clock name.
    key::Vector{T}
    free::Set{Int}
    key_space::F
    prefix::P
    function KeyedRemovalPrefixSearch{T,P}(prefix::P) where {T,P}
        yes = _ -> true
        new{T,P,typeof(yes)}(Dict{T,Int64}(), Vector{T}(), Set{Int}(), yes, prefix)
    end
end


Base.in(kp::KeyedRemovalPrefixSearch, clock_id) = kp.key_space(clock_id)


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
    rng::AbstractRNG, d::Random.SamplerTrivial{KeyedRemovalPrefixSearch{T,P,F}}
    ) where {T,P,F}
    total = sum!(d[])
    choose(d[], rand(rng, Uniform{Float64}(0, total)))
end
