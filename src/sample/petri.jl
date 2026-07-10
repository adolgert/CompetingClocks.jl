"""
    Petri{KeyType,TimeType}(dt=1.0)

A chaos-monkey sampler for smoke-testing. It ignores the clocks'
distributions entirely and picks the next clock to fire uniformly at random
among all currently enabled clocks, advancing time by a fixed step `dt`
(default `1.0`). It's called "Petri" because a Petri net model always chooses
the next event at random.

Because it drives a simulation through improbable event orderings that a
distribution-weighted sampler would rarely produce, it is useful for checking
that a simulation's enable/disable bookkeeping doesn't crash under unlikely
sequences of events. See [`PetriMethod`](@ref) to select it through the
builder interface.
"""
mutable struct Petri{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    time_duration::T
    streams::KeyedStreams{K}
    Petri{K,T}(dt=1.0; seed=_DEFAULT_STREAM_SEED) where {K,T} =
        new(Dict{K,EnablingEntry{K,T}}(), dt, KeyedStreams{K}(seed))
end


# The winner is picked by a single race-level draw (which enabled clock), so
# Petri uses the reserved race stream rather than a per-clock stream.
function clone(propagator::Petri{K,T}) where {K,T}
    c = Petri{K,T}(propagator.time_duration; seed=propagator.streams.seed)
    copy!(c.enabled, propagator.enabled)
    c.streams = copy(propagator.streams)
    return c
end

similar_sampler(propagator::Petri{K,T}) where {K,T} =
    Petri{K,T}(propagator.time_duration; seed=propagator.streams.seed)

function copy_clocks!(dst::Petri, src::Petri)
    copy!(dst.enabled, src.enabled)
    dst.time_duration = src.time_duration
    dst.streams = copy(src.streams)
    return dst
end

rekey_streams!(propagator::Petri, seed) = (rekey_streams!(propagator.streams, seed); propagator)

function next(propagator::Petri{K,T}, when::T) where {K,T}
    kk = keys(propagator.enabled)
    isempty(kk) && return (typemax(T), nothing)
    chosen = rand(race_stream(propagator.streams), kk)
    (when + propagator.time_duration, chosen)
end
