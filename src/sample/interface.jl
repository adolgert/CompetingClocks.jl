using Random: AbstractRNG
using Distributions: UnivariateDistribution
import Base: getindex, keys, length, keytype

export enable!, disable!, next, 
    getindex, keys, length, keytype


"""
    enable!(sampler, clock, distribution, enablingtime, currenttime, RNG)

Tell the sampler to start a clock.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Distributions.UnivariateDistribution`
 * `enablingtime::TimeType` - The zero time for the clock's distribution, in absolute time. Usually equal to `when`.
 * `when::TimeType` - The current time of the simulation.
 * `rng::AbstractRNG` - A random number generator.

"""
function enable!(
    sampler::SSA{K,T},
    clock::K,
    distribution::UnivariateDistribution,
    te::T, # enabling time
    when::T, # current simulation time
    rng::AbstractRNG
    ) where {K,T}
end

"""
    disable!(sampler, clock, when)

Tell the sampler to forget a clock. We include the current simulation time
because some Next Reaction methods use this to optimize sampling.
"""
function disable!(sampler::SSA{K,T}, clock::K, when::T) where {K,T} end

"""
    next(sampler, when, rng)

Ask the sampler for what happens next, in the form of
`(when, which)::Tuple{TimeType,KeyType}`. `rng` is a random number generator.
"""
function next(sampler::SSA{K,T}, when::T, rng::AbstractRNG) where {K,T} end

"""
    getindex(sampler, clock::KeyType)

Return stored state for a particular clock. If the clock does not exist,
a `KeyError` will be thrown.
"""
function Base.getindex(sampler::SSA{K,T}, clock::K) where {K,T} end

"""
    keys(sampler)

Return all stored clocks as a vector.
"""
function Base.keys(sampler::SSA) end

"""
    length(sampler)::Int64

Return the number of stored clocks.
"""
function Base.length(sampler::SSA) end

"""
    keytype(sampler)

Return the type of clock keys.
"""
Base.keytype(::SSA{K,T}) where {K,T} = K
