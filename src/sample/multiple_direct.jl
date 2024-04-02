export MultipleDirect


mutable struct MultipleDirect{SamplerKey,K,Time,Chooser}
    scan::Vector{KeyedPrefixSearch}
    totals::Vector{Time}
    chooser::Chooser
    chosen::Dict{K,Int}
    # Use this indirection from samplerkey to integer so that the
    # list of scan totals is stable, not jumbled by dict keys changing.
    scanmap::Dict{SamplerKey,Int}
end

function MultipleDirect{SamplerKey,K,Time}(
    chooser::Chooser
) where {SamplerKey,K,Time,Chooser <: SamplerChoice{K,SamplerKey}}
    MultipleDirect{SamplerKey,K,Time,Chooser}(
        Vector{KeyedPrefixSearch}(),
        Vector{Time}(),
        chooser,
        Dict{K,SamplerKey}(),
        Dict{SamplerKey,Int}()
    )
end


function reset!(md::MultipleDirect)
    for prefix_search in md.scan
        empty!(prefix_search)
        empty!(md.totals)
        empty!(md.chosen)
        empty!(md.scanmap)
    end
end


function Base.setindex!(
    md::MultipleDirect{SamplerKey,K,Time,Chooser}, keyed_prefix_search, key
    ) where {SamplerKey,K,Time,Chooser}

    push!(md.scan, keyed_prefix_search)
    push!(md.totals, zero(Time))
    md.scanmap[key] = length(md.scan)
end


function enable!(md::MultipleDirect, clock, distribution::Exponential,
    te, when, rng::AbstractRNG)
    if clock ∉ keys(md.chosen)
        which_prefix_search = choose_sampler(md.chooser, clock, distribution)
        scan_idx = md.scanmap[which_prefix_search]
        md.chosen[clock] = scan_idx
    else
        scan_idx = md.chosen[clock]
    end
    keyed_prefix_search = md.scan[scan_idx]
    keyed_prefix_search[clock] = rate(distribution)
end


function disable!(md::MultipleDirect, clock, when)
    which_prefix_search = md.chosen[clock]
    delete!(md.scan[which_prefix_search], clock)
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
function next(md::MultipleDirect, when, rng::AbstractRNG)
    for scan_idx in eachindex(md.scan)
        md.totals[scan_idx] = sum!(md.scan[scan_idx])
    end
    total = sum(md.totals)
    if total > eps(when)
        tau = when + rand(rng, Exponential(1 / total))
        md.totals /= total
        chosen_idx = rand(rng, Categorical(md.totals))
        chosen, hazard_value = rand(rng, md.scan[chosen_idx])
        return (tau, chosen)
    else
        return (typemax(when), nothing)
    end
end
