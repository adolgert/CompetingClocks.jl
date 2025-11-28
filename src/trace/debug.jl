using Base
using Logging
using Random: AbstractRNG
using Distributions: UnivariateDistribution

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
    log::Bool
    DebugWatcher{K,T}(; log=true) where {K,T} =
        new(Vector{EnablingEntry{K,T}}(), Vector{DisablingEntry{K,T}}(), log)
end


clone(ts::DebugWatcher{K,T}) where {K,T} = DebugWatcher{K,T}(; log=ts.log)


function reset!(ts::DebugWatcher)
    ts.log && @debug "reset!"
    empty!(ts.enabled)
    empty!(ts.disabled)
    nothing
end

function copy_clocks!(dst::DebugWatcher{K,T}, src::DebugWatcher{K,T}) where {K,T}
    src.log && @debug "copy!"
    copy!(dst.enabled, src.enabled)
    copy!(dst.disabled, src.disabled)
    return dst
end

function enable!(ts::DebugWatcher{K,T}, clock::K, dist::UnivariateDistribution, te::T, when::T, rng::AbstractRNG) where {K,T}
    ts.log && @debug "enable! $(clock) $(dist) $(te) $(when)"
    push!(ts.enabled, EnablingEntry(clock, dist, te, when))
end


function disable!(ts::DebugWatcher{K,T}, clock::K, when::T) where {K,T}
    ts.log && @debug "disable! $(clock) $(when)"
    push!(ts.disabled, DisablingEntry(clock, when))
end


function fire!(ts::DebugWatcher{K,T}, clock::K, when::T) where {K,T}
    ts.log && @debug "fire! $(clock) $(when)"
    push!(ts.disabled, DisablingEntry(clock, when))
end
