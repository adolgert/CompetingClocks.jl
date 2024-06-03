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


function steploglikelihood(tw::TrackWatcher{K,T}, now, when, which) where {K,T}
    # It's the log of \lambda(t_1) \prod S(t_0, t_1). Here, S is conditional survival, so
    # it's S(t_1)/S(t_0). That's \lambda(t_1) \prod S(t_1)/S(t_0). There isn't a good
    # function to get hazard by itself, so use pdf(t_1)/S(t_1).
    # Take the log, and you get:
    # logpdf(t_1) - logccdf(t_1) + \sum logccdf(t_1) - logccdf(t_0)
    chosen = tw.enabled[which]
    total = logpdf(chosen.distribution, now - chosen.te) - logccdf(chosen.distribution, now - chosen.te)
    for entry in values(tw.enabled)
        total += logccdf(entry.distribution, now - entry.te) - logccdf(entry.distribution, when - entry.te)
    end
    return total
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
