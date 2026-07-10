```@meta
CurrentModule = CompetingClocks
```

# Re-evaluation Couplings

A clock's distribution often changes while the clock is still running: a rate
depends on state, and the state just moved. The clock keeps its age — the
question is what happens to the *random draw* that was in flight. Version 0.4
makes that choice explicit. [`reenable!`](@ref) re-evaluates a
currently-enabled clock's distribution and takes a `coupling` argument
declaring which of two couplings the change implements:

- `:redraw` — discard any retained draw or schedule and draw the remaining
  lifetime fresh, conditioned on the clock's current age.
- `:carry` — deterministic carry: consume *no* randomness; map the retained
  draw through the distribution change by matching conditional survival.

Both couplings produce the same *law* — a statistical test cannot tell them
apart — and they differ exactly where pathwise derivatives care, which is why
the choice belongs to the caller, not the backend. (Calling plain `enable!`
on an already-enabled clock still works, but which coupling it implements is
then chosen silently per backend; `reenable!` exists so nothing is silent.
See the [migration note](migration_04.md) on
this hazard.)

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

Which samplers can carry is a capability:
[`supports_carry`](@ref) is `true` for `CombinedNextReaction` (it retains the
survival uniform outright) and `FirstToFire` (it reconstructs the draw from
its stored firing time), and `false` for `FirstReaction` and the
exponential-only samplers, where `reenable!(..., :carry)` throws an
`ArgumentError` naming the trait. `FirstReaction` accepts either symbol and
treats both as its (moot) redraw-everything re-evaluation, since it redraws
every clock at every `next` anyway.

```julia
using CompetingClocks: CombinedNextReaction, enable!, reenable!
using Distributions

s = CombinedNextReaction{Symbol,Float64}(2024)
enable!(s, :A, Exponential(1 / 0.8), 0.0, 0.0)   # rate 0.8 from time 0

# At time 0.3 the rate becomes 1.5. Carry: no randomness consumed, the
# remaining time rescales by 0.8/1.5, and the schedule moves continuously
# in the new rate.
reenable!(s, :A, Exponential(1 / 1.5), 0.0, 0.3, :carry)
```

## The enabling-time convention: absolute versus relative

`reenable!` carries the same enabling-time flexibility as `enable!`, and the
two API levels anchor it differently. This is the sharpest edge of the verb.

**At the sampler level** the signature is
`reenable!(sampler, clock, dist, te, when, coupling)` with `te` an
**absolute** time, exactly like the low-level `enable!`. Keeping the clock's
age means passing its *original* enabling time:

```julia
reenable!(s, :A, new_dist, 0.0, τ, :carry)   # te = 0.0: the original te, age kept
```

**At the context level** the signature is
`reenable!(ctx, clock, dist, relative_te, coupling)` and `relative_te` is a
**shift from the current time**, exactly like the context-level `enable!`:
the new enabling time is `te = time(ctx) + relative_te`. Keeping the clock's
age therefore means passing `relative_te = original_te − time(ctx)` — a
*nonpositive* number:

```julia
# The clock was enabled at te = 0; the context is now at time τ.
reenable!(ctx, :A, new_dist, -τ, :carry)     # te = τ + (−τ) = 0: age kept
```

The two-argument form `reenable!(ctx, clock, dist, coupling)` uses
`relative_te = 0`, which anchors the distribution at *now* and resets the age
to zero — the re-anchoring case, **not** the age-keeping one. If you mean
"same clock, new rate, keep its age," you must pass the shift explicitly. A
related trap at the model level: under `:carry`, returning the current time
as the enabling time re-anchors the clock and *moves the schedule even when
the distribution is unchanged*. The safe idiom is to hold on to the clock's
first enabling time and always pass it (or the equivalent relative shift)
through the re-evaluation.

## What the observers see

The context-level `reenable!` fans out to watchers identically to `enable!`:
the likelihood watcher closes the old segment (banking its survival) and
opens the new one, and a [`TrajectoryRecorder`](@ref)'s live table switches
to the new `(distribution, te)`, so a subsequently recorded firing
back-calculates its uniform against the *current* segment's distribution.
The coupling choice affects only the sampler's retained draw — the recorded
law and the piecewise likelihood are the same under `:carry` and `:redraw`,
which is the statistical statement that the two couplings agree in law.

## Reference

```@docs
reenable!
```
