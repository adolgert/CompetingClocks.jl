import Base: sum!, push!, setindex!, length
using Random
using Distributions: Uniform


"""
    CumSumPrefixSearch{T}()

This stores hazard rates in order to make it easier for the Direct
method to sample them. This version is the simplest possible, but it can
be faster when there are few hazards enabled. It uses a simple array
and, each time the Direct method samples, this evaluates the cumulative
sum of the array.
"""
mutable struct CumSumPrefixSearch{T<:Real}
    array::Vector{T}
    cumulant::Vector{T}
    dirty::Bool
end


function CumSumPrefixSearch{T}() where {T<:Real}
    CumSumPrefixSearch{T}(Vector{T}(), Vector{T}(), false)
end


function Base.empty!(ps::CumSumPrefixSearch)
    empty!(ps.array)
    empty!(ps.cumulant)
    ps.dirty = false
end

function Base.copy!(dst::CumSumPrefixSearch{T}, src::CumSumPrefixSearch{T}) where {T}
    copy!(dst.array, src.array)
    copy!(dst.cumulant, src.cumulant)
    dst.dirty = src.dirty
    return dst
end

Base.length(ps::CumSumPrefixSearch) = length(ps.array)
time_type(ps::CumSumPrefixSearch{T}) where {T} = T
time_type(ps::Type{CumSumPrefixSearch{T}}) where {T} = T


function Base.push!(ps::CumSumPrefixSearch{T}, value::T) where {T}
    push!(ps.array, value)
    push!(ps.cumulant, zero(T))
    ps.dirty = true
    ps
end


function Base.setindex!(ps::CumSumPrefixSearch{T}, value::T, index) where {T}
    ps.array[index] = value
    ps.dirty = true
    value
end
Base.getindex(ps::CumSumPrefixSearch, index) = ps.array[index]

function Base.sum!(ps::CumSumPrefixSearch{T})::T where {T}
    !ps.dirty && return ps.cumulant[end]
    cumsum!(ps.cumulant, ps.array)
    ps.dirty = false
    ps.cumulant[end]
end


function choose(ps::CumSumPrefixSearch{T}, variate::T) where {T}
    ps.dirty && sum!(ps)
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


function isenabled(ps::CumSumPrefixSearch{T}, clock) where {T}
    if 0 < clock ≤ length(ps.array)
        return ps.array[clock] > zero(T)
    else
        return false
    end
end
