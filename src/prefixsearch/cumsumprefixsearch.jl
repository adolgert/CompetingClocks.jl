import Base: sum!, push!, setindex!, length
using Random
using Distributions: Uniform


"""
    CumSumPrefixSearch{T}()

This stores hazard rates in order to make it easier for the Direct
method to sample them. This version is the dumbest possible, but it can
be faster when there are few hazards enabled. It uses a simple array
and, each time the Direct method samples, this evaluates the cumulative
sum of the array.
"""
struct CumSumPrefixSearch{T<:Real}
	array::Vector{T}
    cumulant::Vector{T}
end


function CumSumPrefixSearch{T}() where {T<:Real}
    CumSumPrefixSearch{T}(Vector{T}(), Vector{T}())
end


function Base.empty!(ps::CumSumPrefixSearch)
    empty!(ps.array)
    empty!(ps.cumulant)
end

function Base.copy!(dst::CumSumPrefixSearch{T}, src::CumSumPrefixSearch{T}) where {T}
    copy!(dst.array, src.array)
    copy!(dst.cumulant, src.cumulant)
end

Base.length(ps::CumSumPrefixSearch) = length(ps.array)
time_type(ps::CumSumPrefixSearch{T}) where {T} = T
time_type(ps::Type{CumSumPrefixSearch{T}}) where {T} = T


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
