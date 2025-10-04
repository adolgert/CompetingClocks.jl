using Base
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
