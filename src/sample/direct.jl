using DataStructures
import Random: rand
import Distributions: Uniform

"""
Classic Direct method for exponential transitions.
This doesn't do any caching of rates.
"""
struct MarkovDirect
end


function next(rm::MarkovDirect, process, when, rng)
    total = 0.0
    cumulative = zeros(Float64, 0)
    keys = Array{Any,1}()
    hazards(process, rng) do clock, distribution, enabled, rng2
        total += params(distribution)[1]
        push!(cumulative, total)
        push!(keys, clock)
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
