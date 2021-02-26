using DataStructures
using Random: rand, AbstractRNG
using Distributions: Uniform, Exponential, params

"""
Classic Direct method for exponential transitions.
This doesn't do any caching of rates.
"""
struct MarkovDirect
end


function next(rm::MarkovDirect, process, when, rng::AbstractRNG)
    total = 0.0
    cumulative = zeros(Float64, 0)
    keys = Array{Any,1}()
    hazards(process, rng) do clock, distribution::Exponential, enabled::Bool
        if enabled
            total += params(distribution)[1]
            push!(cumulative, total)
            push!(keys, clock)
        end
    end

    if total > eps(Float64)
        chosen = searchsortedfirst(cumulative, rand(rng, Uniform(0, total)))
        @assert chosen < length(cumulative) + 1
        return (when - log(rand(rng)) / total, keys[chosen])
    else
        return (Inf, nothing)
    end
end

Observer(fr::MarkovDirect) = (hazard, time, updated, rng) -> nothing


mutable struct DirectCall{T}
    total::Float64
    cumulative::Vector{Float64}
    keys::Vector{T}
end


function DirectCall(::Type{T}) where {T}
    DirectCall{T}(0, zeros(Float64, 0), zeros(T, 0))
end


function zero!(dc::DirectCall{T}) where {T}
    dc.total = 0
    dc.cumulative = zeros(Float64, 0)
    dc.keys = zeros(T, 0)
end


function set_clock!(dc::DirectCall{T}, clock::T, distribution::Exponential, enabled, rng::AbstractRNG) where {T}
    if Base.:(==)(enabled, :Enabled)  # Why is == not found otherwise?
        dc.total += params(distribution)[1]
        push!(dc.cumulative, dc.total)
        push!(dc.keys, clock)
    end  # else it's disabled.
end


function next(dc::DirectCall, when::Float64, rng::AbstractRNG)
    if dc.total > eps(Float64)
        chosen = searchsortedfirst(dc.cumulative, rand(rng, Uniform(0, dc.total)))
        @assert chosen < length(dc.cumulative) + 1
        return (when - log(rand(rng)) / dc.total, dc.keys[chosen])
    else
        return (Inf, nothing)
    end
end
