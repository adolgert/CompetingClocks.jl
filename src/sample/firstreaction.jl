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
