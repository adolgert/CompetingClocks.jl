using Distributions: Exponential, params, truncated
using Random: rand, AbstractRNG

using Logging

export FirstReaction, ChatReaction

"""
    FirstReaction{KeyType,TimeType}()

This is the classic first reaction method for general distributions. 
Every time you sample, this goes to each distribution and asks
when it would fire. Then it takes the soonest and throws out the rest of the
sampled times until the next sample. It can also be very fast when there are only a few clocks to sample.

One interesting property of this sampler is that you can call `next()`
multiple times in order to get a distribution of next firing clocks and their
times to fire.
"""
struct FirstReaction{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    FirstReaction{K,T}() where {K,T<:ContinuousTime} = new(Dict{K,EnablingEntry{K,T}}())
end


clone(fr::FirstReaction{K,T}) where {K,T} = FirstReaction{K,T}()

function _sample_time(entry::EnablingEntry, when::T, rng) where {T}
    dist = entry.te < when ?
           truncated(entry.distribution, when - entry.te, typemax(T)) :
           entry.distribution
    return entry.te + rand(rng, dist)
end


function next(fr::FirstReaction{K,T}, when::T, rng::AbstractRNG) where {K,T}
    isempty(fr.enabled) && return (typemax(T), nothing)
    return minimum(
        ((_sample_time(entry, when, rng), entry.clock) for entry in values(fr.enabled))
    )
end


"""
This sampler can help if it's the first time you're trying a model. It checks
all of the things and uses Julia's logger to communicate them. It samples
using the first reaction algorithm.
"""
mutable struct ChatReaction{K,T}
    # This other class already stores the current set of distributions, so use it.
    enabled::TrackWatcher{K}
    step_cnt::Int64
    enables::Set{K}
    disables::Set{K}
    ChatReaction{K,T}() where {K,T<:ContinuousTime} = new(TrackWatcher{K,T}(), 0, Set{K}(), Set{K}())
end


function enable!(fr::ChatReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T, rng::AbstractRNG) where {K,T}

    # It's OK i fhte clock ∈ keys(fr.enabled.enabled)
    enable!(fr.enabled, clock, distribution, te, when, rng)
    push!(fr.enables, clock)
end


function disable!(fr::ChatReaction{K,T}, clock::K, when::T) where {K,T}
    if clock ∉ fr.enables
        @error "Disabling a clock that was never enabled: $(clock)."
    end
    if clock ∉ keys(fr.enabled.enabled)
        @error "Disabling clock $(clock) that wasn't currently enabled"
    end
    disable!(fr.enabled, clock, when)
    push!(fr.disables, clock)
end


function next(fr::ChatReaction{K,T}, when::T, rng) where {K,T}
    soonest_clock = nothing
    soonest_time = typemax(T)

    for entry::EnablingEntry{K,T} in fr.enabled
        if entry.te < when
            relative_dist = truncated(entry.distribution, when - entry.te, typemax(T))
            putative_time = entry.te + rand(rng, relative_dist)
        else
            putative_time = entry.te + rand(rng, entry.distribution)
        end
        if putative_time < soonest_time
            soonest_clock = entry.clock
            soonest_time = putative_time
        end
    end
    fr.step_cnt += 1
    @debug "Step $(fr.step_cnt) time $(soonest_time) fire $(soonest_clock)"
    return (soonest_time, soonest_clock)
end


Base.haskey(fr::ChatReaction, clock) = haskey(fr.enabled, clock)
