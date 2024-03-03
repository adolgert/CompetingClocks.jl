import Base: sum!, push!, setindex!, length
using Random
using Distributions: Uniform


struct CumSumPrefixSearch{T<:Real}
	array::Vector{T}
    cumulant::Vector{T}
end


function CumSumPrefixSearch{T}() where {T<:Real}
    CumSumPrefixSearch{T}(Vector{T}(), Vector{T}())
end


Base.length(ps::CumSumPrefixSearch) = length(ps.array)


function Base.push!(ps::CumSumPrefixSearch{T}, value::T) where {T}
    push!(ps.array, value)
    push!(ps.cumulant, zero(T))
    ps
end


function Base.setindex!(ps::CumSumPrefixSearch{T}, value::T, index) where {T}
    ps.array[index] = value
    value
end


function Base.sum!(ps::CumSumPrefixSearch{T})::T where {T}
    cumsum!(ps.cumulant, ps.array)
    ps.cumulant[end]
end


function choose(ps::CumSumPrefixSearch{T}, variate::T) where {T}
    index = searchsortedfirst(ps.cumulant, variate)
    return (index, ps.array[index])
end


"""
    rand(rng, sampler::SamplerTrivial{CumSumPrefixSearch})

This method overload allows the machinery of Random to generate random variates
from the CumSumPrefixSearch set of values.
"""
Random.rand(rng::AbstractRNG, d::Random.SamplerTrivial{CumSumPrefixSearch{T}}) where {T} =
    choose(d[], rand(rng, Uniform{T}(zero(T), d[].cumulant[end])))
