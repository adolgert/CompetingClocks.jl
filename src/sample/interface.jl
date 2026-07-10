using Random: AbstractRNG
using Distributions: UnivariateDistribution
import Base: getindex, keys, length, keytype, haskey

"""
    SSA{KeyType,TimeType}

This abstract type represents a stochastic simulation algorithm (SSA). It is
parametrized by the clock ID, or key, and the type used for the time, which
is typically a Float64. The type of the key can be anything you would use
as a dictionary key. This excludes mutable values but includes a wide range
of identifiers useful for simulation. For instance, it could be a `String`,
but it could be a `Tuple{Int64,Int64,Int64}`, so that it indexes into a
complicated simulation state.
"""
abstract type SSA{Key,Time} end


"""
    enable!(sampler, clock, distribution, enablingtime, currenttime)

Tell the sampler to start a clock.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Distributions.UnivariateDistribution`
 * `enablingtime::TimeType` - The zero time for the clock's distribution, in absolute time. Usually equal to `when`.
 * `when::TimeType` - The current time of the simulation.

No random number generator is passed: the sampler OWNS its randomness as a
family of per-clock keyed streams (see [`CompetingClocks.KeyedStreams`](@ref)),
so a draw belongs to `clock`, not to a position in a global call sequence. Two
samplers built from the same seed therefore draw identically for the same clock
regardless of event order.

These times are **absolute** since the start of the simulation. The current time
should be `when`. If you want to shift the distribution so that this event cannot
happen for a little while then choose `enablingtime > when`. If you want to modify
the distribution by shifting it left, then choose `enablingtime < when`. Usually,
`enablingtime == when`. The `truncated()` function from Distributions.jl can make
the left-shift conditioning explicit, but pair it with the shifted enabling time —
`enable!(sampler, clock, truncated(dist; lower=age), when - age, when)` —
because a truncated distribution with `enablingtime == when` measures the
truncated sample from now, pushing every firing at least `age` into the future.

Calling `enable!` on a clock that is ALREADY enabled re-evaluates its
distribution in place and keeps the clock's age. Which pathwise coupling that
implements is chosen SILENTLY by the backend: `CombinedNextReaction` reuses the
retained survival (deterministic carry), while `FirstToFire`, `FirstReaction`,
and the exponential-only samplers redraw at the change. Both agree in law and
differ only where pathwise (IPA) derivatives care. To make the choice EXPLICIT,
call [`reenable!`](@ref)` (..., :carry | :redraw)` instead; this `enable!`
re-enable path is retained for backward compatibility.
"""
function enable!(
    sampler::SSA{K,T},
    clock::K,
    distribution::UnivariateDistribution,
    te::T, # enabling time
    when::T, # current simulation time
) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end

"""
    reenable!(sampler, clock, distribution, te, when, coupling::Symbol)

Re-evaluate the distribution of a clock that is CURRENTLY ENABLED, keeping its
age, and declare EXPLICITLY which pathwise coupling the change implements. This
is the explicit form of what a plain [`enable!`](@ref) on an already-enabled key
does implicitly — plain `enable!` picks the coupling silently per backend, while
`reenable!` makes it a per-call choice so a pathwise (IPA) estimator can insist
on carry and a redraw-at-change model can insist on redraw.

`clock` must already be enabled in `sampler`. Its distribution becomes
`distribution`; `te` is the (possibly shifted) enabling time of the new segment,
carrying the same flexibility as `enable!` — pass the original enabling time to
keep the clock's age, or a shifted one for re-anchoring. `coupling` is one of:

  * `:redraw` — discard any retained draw/schedule and draw the remaining
    lifetime fresh, conditioned on the clock's current age. This is the
    redraw-at-change coupling; it is exactly what a plain `enable!` on an
    already-enabled key does in `FirstToFire` and the exponential-only samplers.
  * `:carry` — deterministic carry: consume NO randomness; map the retained draw
    through the distribution change by matching conditional survival, so the new
    firing age `a_f'` solves `H_new(a_f') = H_new(a) + H_old(a_f) − H_old(a)`
    with `H = −log S` the integrated hazard, `a` the age at the change, and
    `a_f` the old firing age. Only samplers with
    [`supports_carry`](@ref)` == true` implement it; the others throw an
    `ArgumentError`.

This generic method serves every sampler whose plain `enable!`-on-enabled-key IS
the redraw coupling (`FirstToFire`, the exponential-only backends, `Petri`): it
forwards `:redraw` to `enable!` and rejects `:carry` unless the sampler opts in
with `supports_carry`. `CombinedNextReaction` (whose historical `enable!`
behavior IS carry) and `FirstToFire` (which adds a carry from its stored firing
time) override it.
"""
function reenable!(
    sampler::SSA{K,T},
    clock::K,
    distribution::UnivariateDistribution,
    te::T,
    when::T,
    coupling::Symbol,
) where {K,T}
    if coupling === :redraw
        # Plain enable! on an already-enabled key is the redraw-at-change
        # coupling for every sampler that keeps no retained draw. Backends whose
        # enable! instead carries (CombinedNextReaction) override reenable!.
        enable!(sampler, clock, distribution, te, when)
    elseif coupling === :carry
        supports_carry(sampler) || throw(ArgumentError(
            "reenable!(..., :carry) is not supported by $(typeof(sampler)); " *
            "supports_carry(::$(typeof(sampler))) is false. This sampler keeps " *
            "no in-flight draw to carry through a distribution change; use " *
            ":redraw, or a sampler with supports_carry == true."))
        # A carry-supporting sampler must override reenable!; reaching here means
        # the trait was set true without a method.
        error("supports_carry($(typeof(sampler))) is true but reenable! has no " *
              ":carry method for it — this is a package bug.")
    else
        throw(ArgumentError(
            "coupling must be :carry or :redraw, got :$coupling."))
    end
end

"""
    reset!(sampler)

After a sampler is used for a simulation run, it has internal state. This
function resets that internal state to the initial value in preparation
for another sample run.
"""
function reset!(sampler::SSA{K,T}) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    copy_clocks!(destination_sampler, source_sampler)

This copies the state of the source sampler to the destination sampler, replacing
the current state of the destination sampler. This is useful for splitting
techniques where you make copies of a simulation and restart it with different
random number generators.
"""
function copy_clocks!(dst::SSA, src::SSA)
    error("Not implemented for $(typeof(dst))")
end


"""
    clone(sampler)

Make a FULL-STATE deep copy of `sampler`: its clock state, distributions,
schedules/heaps, AND its keyed random streams including every live generator's
state and occurrence counts. Running the clone and the original forward produces
IDENTICAL firing sequences (they are coupled) until one is `rekey_streams!`-ed.
This is the primitive behind clone coupling and common random numbers.

For an EMPTY same-type sampler (an initialized-but-unused copy, e.g. to fill an
array of samplers), use [`similar_sampler`](@ref) instead.
"""
function clone(sampler::SSA{K,T}) where {K,T}
    error("Clone not implemented for $(typeof(sampler))")
end


"""
    similar_sampler(sampler)

Given an existing sampler, make a copy that has the same type and same
constructor options but has no clock data in it. Use this to initialize an array
of samplers or to start a fresh run. Its streams are fresh (empty of live
generators); reseed with [`rekey_streams!`](@ref) if you need a different seed.
This is the old empty-copy meaning of `clone`; `clone` now makes a full-state
coupled copy.
"""
function similar_sampler(sampler::SSA{K,T}) where {K,T}
    error("similar_sampler not implemented for $(typeof(sampler))")
end


"""
    rekey_streams!(sampler, seed)

Re-seed the sampler's keyed random streams to `seed`, forgetting every live
per-key generator and count. After this the sampler draws a fresh, independent
trajectory. Pairing `clone` (a coupled copy) with `rekey_streams!` decouples the
copy; `split!` uses it to give each split copy its own randomness.
"""
function rekey_streams!(sampler::SSA{K,T}, seed) where {K,T}
    error("rekey_streams! not implemented for $(typeof(sampler))")
end


"""
    jitter!(sampler, when)

Takes a sampler that is at a statistically-valid stopping time (this is a technical
term) and resamples all clocks currently-enabled so that this copy of the sampler
will yield different results from another copy of the sampler despite its internal
cache of clock data. The resample draws from the sampler's OWN keyed streams, so
pairing `jitter!` with a preceding [`rekey_streams!`](@ref) is what makes the
resample differ from the original (jittering off identical streams would
reproduce the original draws).
"""
function jitter!(sampler::SSA{K,T}, when::T) where {K,T}
    error("jitter! not implemented for $(typeof(sampler))")
end


"""
    disable!(sampler, clock, when)

Tell the sampler to forget a clock. We include the current simulation time
because some Next Reaction methods use this to optimize sampling.
"""
function disable!(sampler::SSA{K,T}, clock::K, when::T) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end

"""
    fire!(sampler, clock, when)

Tell the sampler that `clock` fired at time `when`.

Firing and disabling are different operations on a clock's draw. Disabling
*censors* the draw: the simulation learns only that the clock had not yet
fired, so a sampler may retain the draw's remaining randomness and reuse it
if the clock is enabled again. Firing *realizes* the draw: the firing time
is observed exactly, no residual randomness remains, and a re-enabled clock
must draw fresh.

The two operations act identically on stored state whenever a sampler
retains no residual draw randomness across state changes. That holds for
every sampler in this package except [`CombinedNextReaction`](@ref), whose
per-clock survival is retained and reused across enabling episodes
(the Anderson/Gibson-Bruck property). This fallback therefore forwards to
[`disable!`](@ref). Override `fire!` for any sampler that retains residual
draw randomness, so that firing consumes the draw completely.
"""
fire!(sampler::SSA{K,T}, clock::K, when::T) where {K,T} =
    disable!(sampler, clock, when)


"""
    next(sampler, when)

Ask the sampler for what happens next, in the form of
`(when, which)::Tuple{TimeType,KeyType}`. The sampler draws from its own keyed
streams; no random number generator is passed.

The returned `(when, which)` is a *reservation*, not a commitment. It is valid
only until the next `enable!`, `disable!`, or `fire!` call changes the sampler's
state. Calling `next` twice without an intervening state change returns valid
reservations, but they are not guaranteed to be identical across samplers:
`FirstReaction` redraws on each call, while `CombinedNextReaction` returns the
same cached reservation. Act on the most recent call's result. The two supported
responses are to fire it—call `fire!` with the returned clock and time—or to
decline it and stop the simulation (the fixed-horizon pattern). There is
deliberately no `peek`.

The `when` argument must not decrease from one call to the next, and must never
advance past a pending firing time without that event being fired. The time of
the most recently fired event is always a safe choice. Within this invariant,
re-querying `next` at a later `when` is legal — this is how `MultiSampler`
composes sub-samplers: the global minimum event time cannot exceed any
sub-sampler's pending minimum, so losing sub-samplers are re-queried
in-contract at the advanced global time. Outside the invariant, samplers
legitimately disagree: `FirstReaction` re-conditions every clock on survival to
`when`, while `CombinedNextReaction` returns putative times cached at enabling.
Within it they agree, because every surviving clock's putative time is at least
`when`.
"""
function next(sampler::SSA{K,T}, when::T) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    getindex(sampler, clock::KeyType)

Return stored state for a particular clock. If the clock does not exist,
a `KeyError` will be thrown. Different samplers have **different stored state.**
"""
function Base.getindex(sampler::SSA{K,T}, clock::K) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    keys(sampler)

Return all stored clocks as a vector.
"""
function Base.keys(sampler::SSA)
    error("Not implemented for $(typeof(sampler))")
end


"""
    length(sampler)::Int64

Return the number of stored clocks.
"""
function Base.length(sampler::SSA)
    error("Not implemented for $(typeof(sampler))")
end


"""
    keytype(sampler)

Return the type of clock keys.
"""
Base.keytype(::SSA{K,T}) where {K,T} = K


"""
    timetype(sampler)

Return the type of clock times.
"""
timetype(::SSA{K,T}) where {K,T} = T


"""
    haskey(sampler, key)

Return a boolean.
"""
function Base.haskey(sampler::SSA{K,T}, key) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    enabled(sampler)

Returns a read-only set of currently-enabled clocks.
"""
function enabled(sampler::SSA{K,T}) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end
