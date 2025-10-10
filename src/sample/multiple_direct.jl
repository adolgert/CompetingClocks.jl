export MultipleDirect


# MultipleDirect is a sampler that contains multiple prefix-search
# data structures.
# Here SamplerKey is an identifier for the sampler.
# K is the type of the Clock key.
# Time is a Float64 or other clock time.
# Chooser is a function that selects a prefix-search given a clock key.
mutable struct MultipleDirect{SamplerKey,K,Time,Chooser}
    scan::Vector{KeyedPrefixSearch} # List of prefix-search data structures.
    totals::Vector{Time}
    chooser::Chooser # Function that selects a prefix-search given a key.
    chosen::Dict{K,Int} # Map from clock key to index of prefix-search in `scan`.
    # Use this indirection from samplerkey to integer so that the
    # list of scan totals is stable, not jumbled by dict keys changing.
    # Map from the identifier for the sampler to the vector of samplers.
    scanmap::Dict{SamplerKey,Int}
end

function MultipleDirect{SamplerKey,K,Time}(
    chooser::Chooser
) where {SamplerKey,K,Time,Chooser<:SamplerChoice{K,SamplerKey}}
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


function Base.copy!(
    dst::MultipleDirect{SamplerKey,K,Time,Chooser},
    src::MultipleDirect{SamplerKey,K,Time,Chooser}
) where {SamplerKey,K,Time,Chooser}
    copy!(dst.scan, src.scan)
    copy!(dst.totals, src.totals)
    dst.chooser = deepcopy(src.chooser)
    copy!(dst.chosen, src.chosen)
    copy!(dst.scanmap, src.scanmap)
    dst
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


Base.haskey(md::MultipleDirect, clock) = false

function Base.haskey(
    md::MultipleDirect{SamplerKey,K,Time,Chooser},
    clock::K
) where {SamplerKey,K,Time,Chooser}
    if haskey(md.chosen, clock)
        return isenabled(md.scan[md.chosen[clock]])
    else
        return false
    end
end


struct MDEnable{C,P} <: AbstractSet{C}
    subset::Vector{P}  # Abstract type for different contained DirectCall.
end

function enabled(md::MultipleDirect{SK,K,T,Choose}) where {SK,K,T,Choose}
    subset = collect(enabled(prefix) for prefix in md.scan)
    MDEnable{K,eltype(subset)}(subset)
end

function Base.iterate(mde::MDEnable)
    isempty(mde.subset) && return nothing
    sub_idx = 1
    res = iterate(mde.subset[sub_idx])
    # This chains iterables from sets that may be empty.
    while res === nothing
        if sub_idx >= lastindex(mde.subset)
            return nothing
        else
            sub_idx += 1
            res = iterate(mde.subset[sub_idx])
        end
    end
    (item, substate) = res
    return (item, (substate, sub_idx))
end


function Base.iterate(mde::MDEnable, (substate, sub_idx))
    res = iterate(mde.subset[sub_idx], substate)
    while res === nothing
        if sub_idx >= lastindex(mde.subset)
            return nothing
        else
            sub_idx += 1
            res = iterate(mde.subset[sub_idx])
        end
    end
    return (res[1], (res[2], sub_idx))
end

# The length() method on a prefix sum is the storage length not the number of enabled clocks.
Base.length(mde::MDEnable) = sum(x -> length(x), mde.subset)
Base.in(x, mde::MDEnable) = any(in(x, subset) for subset in mde.subset)
Base.eltype(::Type{MDEnable{C}}) where {C} = C
isenabled(mde::MultipleDirect, x) = any(isenabled(scan, x) for scan in mde.scan)
