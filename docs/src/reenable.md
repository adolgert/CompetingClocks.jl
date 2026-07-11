```@meta
CurrentModule = CompetingClocks
```

# Re-evaluation and the Sampler's Coupling

A clock's distribution often changes while the clock is still running: a rate
depends on state, and the state just moved. The clock keeps its age — the law
of its remaining firing time is then fully specified by the new distribution,
the enabling time `te`, and the current time. [`reenable!`](@ref) says exactly
that and nothing more: *re-evaluate this clock's distribution to `(dist, te)`.*

```julia
reenable!(sampler, clock, dist, te, when)   # sampler level, te absolute
reenable!(ctx, clock, dist, relative_te)    # context level, shift from now
```

## The coupling is a sampler property

*How* a sampler realizes the change pathwise — what happens to the random draw
that was in flight — is not part of the `reenable!` call. It is a property of
the sampler, fixed at construction, because a choice that leaves the law
unchanged is not a statement about the transition's distribution; it is a
property of how the sampler generates its numbers. It belongs neither to the
model nor to the call site. The two couplings are:

- `:carry` — deterministic carry: consume *no* randomness; map the retained
  draw through the distribution change by matching conditional survival.
- `:redraw` — discard any retained draw or schedule and draw the remaining
  lifetime fresh, conditioned on the clock's current age.

Both couplings produce the same *law* — a statistical test cannot tell them
apart — and they differ exactly where pathwise derivatives and
common-random-number coupling care. The carry-capable samplers,
[`CombinedNextReaction`](@ref) and [`FirstToFire`](@ref), take a `coupling`
keyword at construction, defaulting to `:carry`:

```julia
s = CombinedNextReaction{Symbol,Float64}(2024)                    # carries
s = CombinedNextReaction{Symbol,Float64}(2024; coupling=:redraw)  # redraws
```

The `:carry` default makes these samplers' historical behavior on re-enabling
an enabled clock the explicit default: an existing simulation that never
mentioned a coupling gets the same trajectory it always got. The builder-level
sampler specs forward the same keyword — `NextReactionMethod(coupling=:redraw)`
and `FirstToFireMethod(coupling=:redraw)` — and every other sampler retains no
in-flight draw, so redraw is its only possible behavior: those samplers store
no field and reject `coupling=:carry` at construction with an error naming
[`supports_carry`](@ref). Read a sampler's (or a context's) coupling back with
the accessor [`coupling`](@ref), e.g. `CompetingClocks.coupling(s)`.

## When each coupling is right

**`:carry` is the IPA-safe choice.** Under carry, the new firing time is a
deterministic, continuous function of the new distribution's parameters: a
small change in a rate moves the schedule by a proportionally small amount.
That continuity is precisely what pathwise derivative estimation — IPA,
infinitesimal perturbation analysis, which differentiates an outcome along a
fixed realization of the randomness — requires. A redraw at the change point
would replace the in-flight draw with a fresh one, making the firing time
jump discontinuously in the parameter and destroying the pathwise derivative.
Concretely, carry solves the integrated-hazard matching equation

```
H_new(a_f') = H_new(a) + H_old(a_f) − H_old(a)
```

for the new firing age `a_f'`, where `H = −log S` is the integrated hazard,
`a` is the age at the change, and `a_f` the old firing age. For exponential
old and new distributions this reduces exactly to Gibson–Bruck rate-ratio
scaling of the remaining time, `remaining_new = (rate_old / rate_new) ×
remaining_old`. With an *unchanged* distribution, carry is the identity: the
schedule does not move by a single bit.

**`:redraw` is the modeling choice for a genuinely new mechanism.** If the
state change replaced the physical process counting down — not merely
adjusted its rate — a fresh draw of the remaining lifetime, conditioned on
the age already survived, is the honest model. Redraw consumes one draw from
the clock's own keyed stream, so even a redraw stays reproducible and coupled
across a cloned pair (see [Randomness Ownership](randomness.md)).

## Where carry comes from

The carry map is the coupling induced by Kurtz's random-time-change
representation of a continuous-time jump process: each clock is driven by its
own unit-rate Poisson process, with the clock's internal time given by its
integrated hazard. Holding the Poisson draw fixed while the hazard changes is
exactly the conditional-survival matching above. Anderson (2007) made this
representation explicit for the Next Reaction Method — his Modified Next
Reaction Method, which this package's `CombinedNextReaction` implements in log
(integrated-hazard) space — citing Kurtz for the representation.

Gibson & Bruck (2000) originated the draw-reuse idea for the Next Reaction
Method: they prove reusing a clock's draw across state changes legitimate and
defend it on efficiency grounds (one random number per event). Their machinery,
however, lives on exponential clocks — constant rates (their Theorem 1, the
rate-ratio scaling above) or time-varying rates (their §4.1, an inhomogeneous
Poisson clock). Their Theorem 2 states the general survival-matching transform,
but the paper never exercises a carry across a mid-flight change of a genuinely
non-exponential renewal clock; their §4.2 handles non-exponential waiting times
by decomposing them into exponential sub-steps or drawing the composite time
fresh. The general, correct account for the non-exponential carry this
sampler performs is the Anderson–Kurtz random-time-change one. (Gibson & Bruck
do not claim a conditional redraw biases the law; an *unconditioned* redraw
would.)

## Which samplers can carry

[`supports_carry`](@ref) is `true` for `CombinedNextReaction` (it retains the
survival uniform outright) and `FirstToFire` (it reconstructs the draw from
its stored firing time), and `false` for `FirstReaction` and the
exponential-only samplers, whose constructors reject `coupling=:carry`.
`FirstReaction` redraws every clock at every `next` anyway, so its
`reenable!` simply replaces the stored distribution — the coupling question is
moot for it.

```julia
using CompetingClocks: CombinedNextReaction, enable!, reenable!
using Distributions

s = CombinedNextReaction{Symbol,Float64}(2024)   # coupling=:carry by default
enable!(s, :A, Exponential(1 / 0.8), 0.0, 0.0)   # rate 0.8 from time 0

# At time 0.3 the rate becomes 1.5. The sampler carries: no randomness is
# consumed, the remaining time rescales by 0.8/1.5, and the schedule moves
# continuously in the new rate.
reenable!(s, :A, Exponential(1 / 1.5), 0.0, 0.3)
```

## The enabling-time convention: absolute versus relative

`reenable!` carries the same enabling-time flexibility as `enable!`, and the
two API levels anchor it differently. This is the sharpest edge of the verb.

**At the sampler level** the signature is
`reenable!(sampler, clock, dist, te, when)` with `te` an
**absolute** time, exactly like the low-level `enable!`. Keeping the clock's
age means passing its *original* enabling time:

```julia
reenable!(s, :A, new_dist, 0.0, τ)   # te = 0.0: the original te, age kept
```

**At the context level** the signature is
`reenable!(ctx, clock, dist, relative_te)` and `relative_te` is a
**shift from the current time**, exactly like the context-level `enable!`:
the new enabling time is `te = time(ctx) + relative_te`. Keeping the clock's
age therefore means passing `relative_te = original_te − time(ctx)` — a
*nonpositive* number:

```julia
# The clock was enabled at te = 0; the context is now at time τ.
reenable!(ctx, :A, new_dist, -τ)     # te = τ + (−τ) = 0: age kept
```

The three-argument form `reenable!(ctx, clock, dist)` uses
`relative_te = 0`, which anchors the distribution at *now* and resets the age
to zero — the re-anchoring case, **not** the age-keeping one. If you mean
"same clock, new rate, keep its age," you must pass the shift explicitly. A
related trap at the model level: on a carrying sampler, returning the current
time as the enabling time re-anchors the clock and *moves the schedule even
when the distribution is unchanged*. The safe idiom is to hold on to the
clock's first enabling time and always pass it (or the equivalent relative
shift) through the re-evaluation.

## What the observers see

The context-level `reenable!` fans out to watchers identically to `enable!`:
the likelihood watcher closes the old segment (banking its survival) and
opens the new one, and a [`TrajectoryRecorder`](@ref)'s live table switches
to the new `(distribution, te)`, so a subsequently recorded firing
back-calculates its uniform against the *current* segment's distribution.
The sampler's coupling affects only its retained draw — the recorded law and
the piecewise likelihood are the same under `:carry` and `:redraw`, which is
the statistical statement that the two couplings agree in law.

## Reference

```@docs
reenable!
coupling
validate_coupling
```
