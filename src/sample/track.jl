using Distributions: UnivariateDistribution
export TrackWatcher, DebugWatcher, enable!, disable!, steploglikelihood
export trajectoryloglikelihood, fire!

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


"""
    steploglikelihood(tw::TrackWatcher, now, when, which_fires)

Calculate the log-likelihood of a single step in which the `which_fires`
transition fires next. `now` is the current time. `when` is the time when
`which_fires` happens so when > now. You have to call this before the transition fires so that
it is before transitions are enabled and disabled from the previous step.
"""
function steploglikelihood(tw::TrackWatcher{K,T}, now, when, which_fires) where {K,T}
    # Look for a description of this in docs/notes/distributions.pdf, under log-likelihood.
    if which_fires !== nothing
        chosen = tw.enabled[which_fires]
        if now >= chosen.te
            total = logpdf(chosen.distribution, when - chosen.te)
            if chosen.te < now
                # This time-shifts the pdf, usually seen as f(t,t0) = f(t)/(1-F(t0))
                total -= logccdf(chosen.distribution, now - chosen.te)
            end
        else
            # If a transition fires before it's enabled, that's impossible.
            total = -NaN
        end
    else
        total = zero(Float64)
    end
    for (key, entry) in pairs(tw.enabled)
        if key !== which_fires
            if when > entry.te
                total += logccdf(entry.distribution, when - entry.te)
                if now > entry.te
                    total -= logccdf(entry.distribution, now - entry.te)
                end
            end
        end
    end
    return total
end


mutable struct TrajectoryWatcher{K,T}
    track::TrackWatcher{K,T}
    loglikelihood::Float64
    curtime::Float64
    TrajectoryWatcher{K,T}() where {K,T} = new(TrackWatcher{K,T}(), zero(Float64), zero(Float64))
end
export TrajectoryWatcher


function trajectoryloglikelihood(tw::TrajectoryWatcher)
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    when = tw.curtime
    remaining = zero(Float64)
    for entry in values(tw.track.enabled)
        if when > entry.te
            remaining += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                remaining -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
    end

    return tw.loglikelihood + remaining
end


reset!(tw::TrajectoryWatcher) = (reset!(tw.track); tw.loglikelihood=zero(Float64); nothing)
function Base.copy!(dst::TrajectoryWatcher{K,T}, src::TrajectoryWatcher{K,T}) where {K,T}
    copy!(dst.track, src.track)
    dst.loglikelihood = src.loglikelihood
end

Base.iterate(ts::TrajectoryWatcher) = iterate(ts.track)
Base.iterate(ts::TrajectoryWatcher, i::Int64) = iterate(ts.track, i)
Base.length(ts::TrajectoryWatcher) = length(ts.track)


function enable!(ts::TrajectoryWatcher{K,T}, clock::K, dist::UnivariateDistribution, te::T, when::T, rng) where {K,T}
    enable!(ts.track, clock, dist, te, when, rng)
end


function disable!(ts::TrajectoryWatcher{K,T}, clock::K, when) where {K,T}
    entry = get(ts.track.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
        disable!(ts.track, clock, when)
    end
end


function fire!(ts::TrajectoryWatcher{K,T}, clock::K, when) where {K,T}
    entry = get(ts.track.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logpdf(entry.distribution, when - entry.te)
        end
        # Adjust for an enabling time that was shifted left.
        if entry.when > entry.te
            ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
        end
        disable!(ts.track, clock, when)
    end
    ts.curtime = when
end


function steploglikelihood(tw::TrajectoryWatcher{K,T}, now, when, which_fires) where {K,T}
    steploglikelihood(tw.track, now, when, which_fires)
end


"""
    DebugWatcher()

For debugging, it helps to have visibility into the simulation. This Watcher
records everything that is enabled or disabled as a list of all enables and
all disabled. It's the complete event history, and you can think of it as
the filtration for the process going forward.

```
watcher = DebugWatcher{String}()
# enable and disable some things.
(watcher.enabled[1].clock,
watcher.enabled[1].distribution,
watcher.enabled[1].te,
watcher.enabled[1].when)
```
"""
mutable struct DebugWatcher{K,T}
    enabled::Vector{EnablingEntry{K,T}}
    disabled::Vector{DisablingEntry{K,T}}
    DebugWatcher{K,T}() where {K,T}=new(Vector{EnablingEntry{K,T}}(), Vector{DisablingEntry{K,T}}())
end


reset!(ts::DebugWatcher) = (empty!(ts.enabled); empty!(ts.disabled); nothing)

function Base.copy!(dst::DebugWatcher{K,T}, src::DebugWatcher{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    copy!(dst.disabled, src.disabled)
end

function enable!(ts::DebugWatcher{K,T}, clock::K, dist::UnivariateDistribution, te, when, rng) where {K,T}
    push!(ts.enabled, EnablingEntry(clock, dist, te, when))
end


function disable!(ts::DebugWatcher{K,T}, clock::K, when) where {K,T}
    push!(ts.disabled, DisablingEntry(clock, when))
end
