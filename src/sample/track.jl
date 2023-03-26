using Distributions: UnivariateDistribution
export TrackWatcher, DebugWatcher, enable!, disable!

# A Watcher has an enable!() and a disable!() function but lacks
# the next() function that a Sampler has. You can attach a watcher
# to a model in order to provide more information about active
# clocks.

struct EnablingEntry{T}
    clock::T
    distribution::UnivariateDistribution
    te::Float64
    when::Float64
end


struct DisablingEntry{T}
    clock::T
    when::Float64
end


"""
    TrackWatcher()

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
mutable struct TrackWatcher{T}
    enabled::Dict{T,EnablingEntry{T}}
    TrackWatcher{T}() where {T}=new(Dict{T,EnablingEntry{T}}())
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


function enable!(ts::TrackWatcher{T}, clock::T, dist::UnivariateDistribution, te, when, rng) where {T}
    ts.enabled[clock] = EnablingEntry{T}(clock, dist, te, when)
end


function disable!(ts::TrackWatcher{T}, clock::T, when) where {T}
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
mutable struct DebugWatcher{T}
    enabled::Vector{EnablingEntry{T}}
    disabled::Vector{DisablingEntry{T}}
    DebugWatcher{T}() where {T}=new(Vector{EnablingEntry{T}}(), Vector{DisablingEntry{T}}())
end


function enable!(ts::DebugWatcher{T}, clock::T, dist::UnivariateDistribution, te, when, rng) where {T}
    push!(ts.enabled, EnablingEntry(clock, dist, te, when))
end


function disable!(ts::DebugWatcher{T}, clock::T, when) where {T}
    push!(ts.disabled, DisablingEntry(clock, when))
end
