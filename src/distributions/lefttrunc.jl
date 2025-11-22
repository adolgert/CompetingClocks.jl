# This file adds to definitions in the Distributions.Truncated struct.
# It specializes truncations for particular distributions in order to
# improve speed and accuracy.
using Distributions
using Random
using Base: rand


"""
Drawing a random number from a left-truncated exponential is particularly simple.
"""
function Base.rand(rng::AbstractRNG, d::Truncated{<:Exponential,Continuous})
    @assert isinf(d.upper) || d.upper === nothing
    return d.lower + rand(rng, d.untruncated)
end
