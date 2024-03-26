using Random: AbstractRNG
using Distributions: UnivariateDistribution
import Base: getindex, keys, length, keytype

export enable!, disable!, next, 
    getindex, keys, length, keytype


function enable!(
    sampler::SSA{K,T},
    clock::K,
    distribution::UnivariateDistribution,
    te::T, # enabling time
    when::T, # current simulation time
    rng::AbstractRNG
    ) where {K,T}
end

function disable!(sampler::SSA{K,T}, clock::K, when::T) where {K,T} end

function next(sampler::SSA{K,T}, when::T, rng::AbstractRNG) where {K,T} end

"""
    Return stored state for a particular clock. If the clock does not exist,
a `KeyError` will be thrown.
"""
function Base.getindex(sampler::SSA{K,T}, clock::K) where {K,T} end

"""
    Return all stored clocks as a vector.
"""
function Base.keys(sampler::SSA) end

"""
    Return the number of stored clocks.
"""
function Base.length(sampler::SSA) end

"""
    Return the type of clock keys.
"""
Base.keytype(::SSA{K,T}) where {K,T} = K
