using Random: AbstractRNG
using Distributions: UnivariateDistribution
import Base: getindex, keys, length, keytype

export enable!, disable!, next, 
    getindex, keys, length, keytype


function enable!(
    sampler::SSA,
    clock,
    distribution::UnivariateDistribution,
    te, # enabling time
    when, # current simulation time
    rng::AbstractRNG
    )
end

function disable!(sampler, clock, when::Float64) end

function next(sampler, when::Float64, rng::AbstractRNG) end

"""
    Return the time associated with a clock. If the clock does not exist,
a `KeyError` will be thrown.
"""
function Base.getindex(sampler::S, clock::K) where {K, T, S<:SSA{K,T}} end

"""
    Return all stored clocks as a vector.
"""
function Base.keys(sampler::S) where {S<:SSA} end

"""
    Return the number of stored clocks.
"""
function Base.length(sampler::S) where {S<:SSA} end

"""
    Return the type of clock keys.
"""
Base.keytype(::SSA{K,T}) where {K,T} = K