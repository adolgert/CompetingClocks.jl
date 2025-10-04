using Base
using Distributions: UnivariateDistribution
export TrackWatcher, DebugWatcher, enable!, disable!, steploglikelihood
export trajectoryloglikelihood, fire!, absolute_enabling

# A Watcher has an enable!() and a disable!() function but lacks
# the next() function that a Sampler has. You can attach a watcher
# to a model in order to provide more information about active
# clocks.

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
mutable struct TrackWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    TrackWatcher{K,T}() where {K,T}=new(Dict{K,EnablingEntry{K,T}}())
end

function absolute_enabling(dst::TrackWatcher{K,T}, clock::K) where {K,T}
    return dst.enabled[clock].when
end

reset!(ts::TrackWatcher) = (empty!(ts.enabled); nothing)

function Base.copy!(dst::TrackWatcher{K,T}, src::TrackWatcher{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
end

function Base.iterate(ts::TrackWatcher)
    return iterate(values(ts.enabled))
end

function Base.iterate(ts::TrackWatcher, i::Int64)
    return iterate(values(ts.enabled), i)
end


function Base.length(ts::TrackWatcher)
    return length(ts.enabled)
end


function enable!(ts::TrackWatcher{K,T}, clock::K, dist::UnivariateDistribution, te, when, rng) where {K,T}
    ts.enabled[clock] = EnablingEntry{K,T}(clock, dist, te, when)
end


function disable!(ts::TrackWatcher{K,T}, clock::K, when) where {K,T}
    if haskey(ts.enabled, clock)
        delete!(ts.enabled, clock)
    end
end

isenabled(ts::TrackWatcher{K,T}, clock::K) where {K,T} = haskey(ts.enabled, clock)
isenabled(ts::TrackWatcher{K,T}, clock) where {K,T} = false


"""
    steploglikelihood(tw::TrackWatcher, now, when_fires, which_fires)

Calculate the log-likelihood of a single step in which the `which_fires`
transition fires next. `now` is the current time. `when_fires` is the time when
`which_fires` happens so `when > now`. You have to call this before the transition fires so that
it is before transitions are enabled and disabled from the previous step.
"""
function steploglikelihood(tw::TrackWatcher{K,T}, t0, t, which_fires) where {K,T}
    # Look for a description of this in docs/notes/distributions.pdf, under log-likelihood.
    @assert t >= t0
    total = zero(Float64)
    for (key, entry) in pairs(tw.enabled)
        if key == which_fires
            if t >= entry.te
                total += logpdf(entry.distribution, t - entry.te)
                if t0 > entry.te
                    # This time-shifts the pdf, usually seen as f(t,t0) = f(t)/(1-F(t0))
                    total -= logccdf(entry.distribution, t0 - entry.te)
                end
            else
                # If a transition fires before it's enabled, that's impossible.
                total = -Inf
            end
        else
            if t > entry.te
                total += logccdf(entry.distribution, t - entry.te)
                if t0 > entry.te
                    total -= logccdf(entry.distribution, t0 - entry.te)
                end
            end
        end
    end
    return total
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

function enable!(propagator::MemorySampler, clock, distribution, te, when, rng)
    enable!(propagator.track, clock, distribution, te, when, rng)
    enable!(propagator.sampler, clock, distribution, te, when, rng)
end

function disable!(propagator::MemorySampler, clock, when)
    disable!(propagator.track, clock, when)
    disable!(propagator.sampler, clock, when)
end

function Base.getindex(propagator::MemorySampler, clock)
    getindex(propagator.sampler, clock)
end


