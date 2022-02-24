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


struct DirectCall{T}
    key::Dict{T, Int64}
    propensity::Vector{Float64}
end


function DirectCall(::Type{T}) where {T}
    DirectCall{T}(Dict{T, Int64}(), zeros(Float64, 0))
end


function set_clock!(dc::DirectCall{T}, clock::T, distribution::Exponential,
        enabled::Symbol, rng::AbstractRNG) where {T}
    if Base.:(==)(enabled, :Enabled)
        hazard = params(distribution)[1]
        if !haskey(dc.key, clock)
            dc.key[clock] = length(push!(dc.propensity, hazard))
        else
            dc.propensity[dc.key[clock]] = hazard
        end
    elseif Base.:(==)(enabled, :Changed)  # Why is == not found otherwise?
        dc.propensity[dc.key[clock]] = params(distribution)[1]
    else  # else it's disabled.
        dc.propensity[dc.key[clock]] = 0.0
    end
end


function next(dc::DirectCall, when::Float64, rng::AbstractRNG)
    total = sum(dc.propensity)
    if total > eps(Float64)
        chosen = searchsortedfirst(cumsum(dc.propensity), rand(rng, Uniform(0, total)))
        @assert chosen < length(dc.propensity) + 1
        key_name = [x for (x, y) in pairs(dc.key) if Base.:(==)(y, chosen)][1]
        return (when - log(rand(rng)) / total, key_name)
    else
        return (Inf, nothing)
    end
end
