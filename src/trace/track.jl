using Base
using Random: AbstractRNG
using Distributions: UnivariateDistribution
export TrackWatcher, DebugWatcher, enable!, disable!, steploglikelihood
export trajectoryloglikelihood, fire!, absolute_enabling

# A Watcher has an enable!() and a disable!() function but lacks
# the next() function that a Sampler has. You can attach a watcher
# to a model in order to provide more information about active
# clocks.
#
# `when`` and `te`` are in absolute times.
# Here are the possible cases:
# Let's call the simulation time `now`
# At the moment of creation `when=now`. `te` can be less than, the same as, or greater than `when`.
# For later calls, `when < now` strictly less than. So you can have
#  * `te<=when < now`,
#  * `when<=te<now`,
#  * `when<te=now`, or
#  * `when<now<te`.
#
struct EnablingEntry{K,T}
    clock::K
    distribution::UnivariateDistribution
    te::T    # The zero-point, in absolute time, for the distribution.
    when::T  # When the distribution was enabled.
end


struct DisablingEntry{K,T}
    clock::K
    when::T
end


abstract type EnabledWatcher{K,T} <: SSA{K,T} end

"""
    TrackWatcher{K,T}()

This Watcher doesn't sample. It records everything enabled.
You can iterate over enabled clocks with a for-loop. If we think of the
model as providing changes in which transitions are enabled or disabled, this
Watcher accumulates those changes to provide a consistent list of all enabled
transitions. Together, a model and this Watcher provide the Semi-Markov core
matrix, or the row of it that is currently known.

```julia
for entry in tracker
    entry.clock
    entry.distribution
    entry.te
    entry.when
end
```
"""
mutable struct TrackWatcher{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    TrackWatcher{K,T}() where {K,T} = new(Dict{K,EnablingEntry{K,T}}())
end

function absolute_enabling(dst::EnabledWatcher{K,T}, clock::K) where {K,T}
    return dst.enabled[clock].when
end

reset!(ts::EnabledWatcher) = (empty!(ts.enabled); nothing)

function copy_clocks!(dst::EnabledWatcher{K,T}, src::EnabledWatcher{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    return dst
end

Base.keys(ts::EnabledWatcher) = keys(ts.enabled)
Base.getindex(ts::EnabledWatcher{K}, key::K) where {K} = getindex(ts.enabled, key)
Base.haskey(ts::EnabledWatcher{K}, key::K) where {K} = haskey(ts.enabled, key)

function Base.iterate(ts::EnabledWatcher)
    return iterate(values(ts.enabled))
end

function Base.iterate(ts::EnabledWatcher, i::Int64)
    return iterate(values(ts.enabled), i)
end


function Base.length(ts::EnabledWatcher)
    return length(ts.enabled)
end


function enable!(ts::EnabledWatcher{K,T}, clock::K, dist::UnivariateDistribution, te::T, when::T, rng::AbstractRNG) where {K,T}
    haskey(ts.enabled, clock) && disable!(ts, clock, when)
    ts.enabled[clock] = EnablingEntry{K,T}(clock, dist, te, when)
end

fire!(ts::EnabledWatcher{K,T}, clock::K, when::T) where {K,T} = disable!(ts, clock, when)

function disable!(ts::EnabledWatcher{K,T}, clock::K, when::T) where {K,T}
    if haskey(ts.enabled, clock)
        delete!(ts.enabled, clock)
    end
end

isenabled(ts::EnabledWatcher, clock) = haskey(ts.enabled, clock)
enabled(ts::EnabledWatcher) = keys(ts.enabled)

function _steploglikelihood(enabled, t0, t, which_fires)
    # Look for a description of this in docs/notes/distributions.pdf, under log-likelihood.
    @assert t >= t0
    return sum(
        function (entry)
            t < entry.te && return (entry.clock == which_fires) ? -Inf : zero(t)
            fired = if entry.clock == which_fires
                logpdf(entry.distribution, t - entry.te)
            else
                logccdf(entry.distribution, t - entry.te)
            end
            base = (t0 > entry.te) ? logccdf(entry.distribution, t0 - entry.te) : zero(t)
            fired - base
        end,
        enabled
    )
end


"""
    steploglikelihood(tw::TrackWatcher, now, when_fires, which_fires)

Calculate the log-likelihood of a single step in which the `which_fires`
transition fires next. `now` is the current time. `when_fires` is the time when
`which_fires` happens so `when > now`. You have to call this before the transition fires so that
it is before transitions are enabled and disabled from the previous step.
"""
function steploglikelihood(tw::EnabledWatcher{K,T}, t0, t, which_fires) where {K,T}
    _steploglikelihood(values(tw.enabled), t0, t, which_fires)
end


mutable struct MemorySampler{S,K,T}
    sampler::S
    track::TrackWatcher{K,T}
    curtime::T
end

function MemorySampler(sampler::Sampler) where {Sampler}
    K = keytype(sampler)
    T = timetype(sampler)
    MemorySampler{Sampler,K,T}(
        sampler, TrackWatcher{K,T}(), zero(T)
    )
end

export MemorySampler

keytype(propagator::MemorySampler{S,K,T}) where {S,K,T} = K

function absolute_enabling(propagator::MemorySampler, clock)
    return absolute_enabling(propagator.track, clock)
end

function next(propagator::MemorySampler, when, rng)
    next(propagator.sampler, when, rng)
end

function enable!(propagator::MemorySampler{S,K,T}, clock::K, distribution::UnivariateDistribution, te::T, when::T, rng::AbstractRNG) where {S,K,T}
    enable!(propagator.track, clock, distribution, te, when, rng)
    enable!(propagator.sampler, clock, distribution, te, when, rng)
end

function disable!(propagator::MemorySampler{S,K,T}, clock::K, when::T) where {S,K,T}
    disable!(propagator.track, clock, when)
    disable!(propagator.sampler, clock, when)
end

function Base.getindex(propagator::MemorySampler, clock)
    getindex(propagator.sampler, clock)
end
