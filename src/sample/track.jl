using Distributions: UnivariateDistribution
export TrackWatcher, DebugWatcher, enable!, disable!

# A Watcher has an enable!() and a disable!() function but lacks
# the next() function that a Sampler has. You can attach a watcher
# to a model in order to provide more information about active
# clocks.

struct EnablingEntry{K,T}
    clock::K
    distribution::UnivariateDistribution
    te::T
    when::T
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


function enable!(ts::DebugWatcher{K,T}, clock::K, dist::UnivariateDistribution, te, when, rng) where {K,T}
    push!(ts.enabled, EnablingEntry(clock, dist, te, when))
end


function disable!(ts::DebugWatcher{K,T}, clock::K, when) where {K,T}
    push!(ts.disabled, DisablingEntry(clock, when))
end
