

struct MultipleDirect{T}
    scan::Vector{KeyedPrefixSearch}
    totals::Vector{Float64}
end


function enable!(dc::MultipleDirect{T}, clock::T, distribution::Exponential,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}
    for owner in dc.scan
        if clock ∈ owner
            owner[clock] = rate(distribution)
            return
        end
    end
    throw(KeyError("Cannot find sampler for clock id $clock"))
end


function disable!(dc::MultipleDirect{T}, clock::T, when::Float64) where {T}
    for owner in dc.scan
        if clock ∈ owner
            delete!(owner, clock)
            return
        end
    end
    throw(KeyError("Cannot find sampler for clock id $clock"))
end


"""
    next(multiple_direct, when, rng)

Selects the next transition to fire and when it fires.

There are two main algorithms for this selection. This implementation handles
the case when there are a lot of clocks or when some clocks have much
smaller hazards. It first draws a random number to choose which subset
of hazards will be used, and then it asks that subset to draw a random
number to choose which particular hazard is used. When there are many hazards,
it is possible that a random number generator will _never_ choose a particular
value because there is no guarantee that a random number generator covers every
combination of bits. Using more draws decreases the likelihood of this problem.
"""
function next(dc::MultipleDirect, when::Float64, rng::AbstractRNG)
    for scan_idx in eachindex(dc.scan)
        dc.totals[scan_idx] = sum!(dc.scan[scan_idx])
    end
    total = sum(dc.totals)
    if total > eps(Float64)
        tau = when + rand(rng, Exponential(1 / total))
        dc.totals /= total
        chosen_idx = rand(rng, dc.totals)
        chosen, hazard_value = rand(rng, dc.scan[chosen_idx])
        return (tau, chosen)
    else
        return (Inf, nothing)
    end
end
