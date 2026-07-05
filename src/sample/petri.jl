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
    Petri{K,T}(dt=1.0) where {K,T} = new(Dict{K,EnablingEntry{K,T}}(), dt)
end


clone(propagator::Petri{K,T}) where {K,T} = Petri{K,T}(propagator.time_duration)


function copy_clocks!(dst::Petri, src::Petri)
    copy!(dst.enabled, src.enabled)
    dst.time_duration = src.time_duration
    return dst
end

function next(propagator::Petri{K,T}, when::T, rng::AbstractRNG) where {K,T}
    kk = keys(propagator.enabled)
    isempty(kk) && return (typemax(T), nothing)
    chosen = rand(rng, kk)
    (when + propagator.time_duration, chosen)
end
