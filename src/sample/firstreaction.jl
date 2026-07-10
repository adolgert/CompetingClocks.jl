using Distributions: Exponential, params, truncated
using Random: rand, AbstractRNG

using Logging

"""
    FirstReaction{KeyType,TimeType}()

This is the classic first reaction method for general distributions. 
Every time you sample, this goes to each distribution and asks
when it would fire. Then it takes the soonest and throws out the rest of the
sampled times until the next sample. It can also be very fast when there are only a few clocks to sample.

One interesting property of this sampler is that you can call `next()`
multiple times in order to get a distribution of next firing clocks and their
times to fire.
"""
mutable struct FirstReaction{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    streams::KeyedStreams{K}
    FirstReaction{K,T}(seed=_DEFAULT_STREAM_SEED) where {K,T<:ContinuousTime} =
        new(Dict{K,EnablingEntry{K,T}}(), KeyedStreams{K}(seed))
end


# A full-state clone: the enabled table AND the keyed streams (generator states,
# counts) are copied, so the clone and original race the same clocks to the same
# times until one is rekeyed. FirstReaction redraws every clock at every `next`,
# so the coupling shows up on the very next call.
function clone(fr::FirstReaction{K,T}) where {K,T}
    c = FirstReaction{K,T}(fr.streams.seed)
    copy!(c.enabled, fr.enabled)
    c.streams = copy(fr.streams)
    return c
end

similar_sampler(fr::FirstReaction{K,T}) where {K,T} = FirstReaction{K,T}(fr.streams.seed)

# copy_clocks! carries both the enabled table and the stream state (the generic
# EnabledWatcher copy_clocks! copies only `enabled`, which would silently couple
# this sampler's future randomness to the source's).
function copy_clocks!(dst::FirstReaction{K,T}, src::FirstReaction{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    dst.streams = copy(src.streams)
    return dst
end

rekey_streams!(fr::FirstReaction, seed) = (rekey_streams!(fr.streams, seed); fr)

# Each clock's conditional-survival redraw comes from THAT clock's own stream, so
# the redraw is per-(clock, occurrence) and independent of event order.
function _sample_time(ks::KeyedStreams, entry::EnablingEntry, when::T) where {T}
    dist = entry.te < when ?
           truncated(entry.distribution, when - entry.te, typemax(T)) :
           entry.distribution
    return entry.te + rand(stream_for!(ks, entry.clock), dist)
end


function next(fr::FirstReaction{K,T}, when::T) where {K,T}
    isempty(fr.enabled) && return (typemax(T), nothing)
    return minimum(
        ((_sample_time(fr.streams, entry, when), entry.clock) for entry in values(fr.enabled))
    )
end


# --- estimator-facing verbs -------------------------------------------------
# FirstReaction is a conditional-redraw backend: it stores only enabling times
# and distributions and redraws every clock's remaining lifetime, conditioned on
# its age, at each `next`. That makes it a natural home for `enabled_ages` (the
# enabling-time table it already keeps IS the answer) and for `force_fire!`
# (forcing is free — a loser's next conditional draw is already the survival-
# conditioned one), but not for `retained_draw`: it keeps no draw to hand back.

supports_enabled_ages(::Type{<:FirstReaction}) = true
supports_force(::Type{<:FirstReaction}) = true

# The enabling-time table, as (key, te) pairs. `entry.clock` is the key.
enabling_times(fr::FirstReaction) =
    ((entry.clock, entry.te) for entry in values(fr.enabled))

"""
    force_fire!(fr::FirstReaction, key, tstar)

Fire `key` at `tstar`. The body is just `fire!` bookkeeping — remove the fired
clock — because every draw this backend makes is conditional on age at the
moment of drawing, so each surviving clock's next `next` already redraws it
conditioned on survival past `tstar`. No random number is consumed.
"""
function force_fire!(fr::FirstReaction{K,T}, key::K, tstar::T) where {K,T}
    haskey(fr.enabled, key) || throw(KeyError(key))
    delete!(fr.enabled, key)
    nothing
end


"""
    reenable!(fr::FirstReaction, clock, distribution, te, when, coupling)

Re-evaluate `clock`'s distribution mid-flight. For FirstReaction the coupling
choice is MOOT: it retains no schedule and redraws every clock's remaining
lifetime, conditioned on age, at every `next`, so there is no in-flight draw to
either carry or discard — this is the redraw-everything coupling. Both `:carry`
and `:redraw` therefore resolve to the same action, replacing the stored
distribution and enabling time, and both are accepted so a driver can issue the
same `reenable!` call against any backend. (Correspondingly `supports_carry` is
`false`: FirstReaction cannot move a firing time CONTINUOUSLY in a parameter,
because it holds no firing time to move — the property carry exists to provide.)
"""
function reenable!(
    fr::FirstReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T, coupling::Symbol) where {K,T}
    haskey(fr.enabled, clock) || throw(ArgumentError(
        "reenable! needs clock $clock currently enabled; use enable! to start a " *
        "clock that is not enabled."))
    (coupling === :carry || coupling === :redraw) || throw(ArgumentError(
        "coupling must be :carry or :redraw, got :$coupling."))
    enable!(fr, clock, distribution, te, when)
    nothing
end
