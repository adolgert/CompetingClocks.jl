```@meta
CurrentModule = CompetingClocks
```

# Estimator-Facing Verbs

Simulating a generalized semi-Markov process (GSMP — a continuous-time model
in which many clocks race and the earliest firing changes the state) needs
only four verbs: `enable!`, `disable!`, `next`, and `fire!`. Estimating the
*derivative* of a simulation's output with respect to model parameters needs
more. A derivative estimator must be able to ask a live sampler which clocks
are enabled and how long each has been running, to read the random draw
behind a scheduled firing, and to impose a firing of its own choosing. This
page describes those verbs, the capability traits that gate them, and which
estimator needs which verb.

## The three verbs

### `enabled_ages`: which clocks are racing, and how old are they

[`enabled_ages(sampler, when)`](@ref enabled_ages) returns every
currently-enabled clock paired with its age `when - te`, where `te` is the
clock's enabling time. The result is **sorted by key**, and that order is
load-bearing: a branching estimator builds a probability-mass function over
the enabled clocks and indexes it by position, so two samplers — or a sampler
and its clone — must present the same clocks in the same order. See
[Contract and Invariants](contract.md).

The same query is available on a non-delayed context as
`enabled_ages(ctx)` (at the context's current time) or
`enabled_ages(ctx, when)`, so an estimator working at the context level need
not reach into the raw sampler.

```julia
using CompetingClocks: FirstToFire, enable!, enabled_ages
using Distributions

s = FirstToFire{Symbol,Float64}(4040)
enable!(s, :z, Weibull(2.0, 1.5), -0.4, 1.0)   # left-shifted: aged before observed
enable!(s, :a, Weibull(1.7, 1.0), 0.0, 1.0)
enable!(s, :m, Exponential(1.2), 0.5, 1.0)

enabled_ages(s, 1.0)
# 3-element Vector{Tuple{Symbol, Float64}}:
#  (:a, 1.0)
#  (:m, 0.5)
#  (:z, 1.4)     ← age exceeds elapsed time; the shift is part of the age
```

An age can exceed the elapsed simulation time: a clock enabled with its
`te` shifted into the past was already old when the sampler first saw it, and
`enabled_ages` reports that true age.

### `retained_draw`: the uniform behind a scheduled firing

[`retained_draw(sampler, key)`](@ref retained_draw) returns a named tuple
`(te, u, logu)`: the clock's enabling time and the survival-space uniform
behind its currently scheduled tentative firing time, satisfying the
retained-draw identity

```
tentative_time == te + invlogccdf(dist, logu),    u == exp(logu)
```

This is the same identity every recorded firing obeys (see
[Recording Trajectories](recording.md)), reported here for a draw that has
*not fired yet*. Only a sampler that retains draw randomness across state
changes can answer, which in this package is `CombinedNextReaction` alone.
Note that the recorded event payload, not this live query, is what the
estimators surveyed so far actually consumed — the pathwise estimator reads
`(te, u)` from the firing record. The live query exists for estimators that
must inspect a schedule mid-race.

### `force_fire!`: impose a firing

[`force_fire!(sampler, key, tstar)`](@ref force_fire!) fires the *chosen*
clock at the *chosen* time, regardless of which clock would have won the
race. This is the branch step of a weak-derivative (measure-valued) estimator:
having computed that the derivative weight requires the trajectory in which
clock `key` fires at `tstar`, the estimator imposes exactly that firing. The
fired clock is consumed exactly as `fire!` consumes it, and every surviving
clock ends up distributed by its lifetime law conditioned on survival past
`tstar`.

The context-level form `force_fire!(ctx, clock, tstar)` additionally advances
the context's time to `tstar` and fans out to the watchers exactly like
`fire!`, so an attached [`TrajectoryRecorder`](@ref) records the forced firing
— with a well-defined uniform `u = ccdf(dist, tstar - te)` — and a replay
cannot tell an imposed firing from a raced one.

```julia
using CompetingClocks: CombinedNextReaction, enable!, force_fire!, next
using Distributions

s = CombinedNextReaction{Symbol,Float64}(70_707)
enable!(s, :x, Weibull(1.7, 1.0), 0.0, 0.0)
enable!(s, :y, Exponential(0.8), 0.0, 0.0)
enable!(s, :z, Weibull(2.0, 1.5), -0.4, 0.0)

force_fire!(s, :y, 0.6)     # :y fires at 0.6, whatever the race said
t2, w2 = next(s, 0.6)       # survivors are conditioned on survival past 0.6
@assert w2 in (:x, :z) && t2 > 0.6
```

Two semantic points matter. First, no random number generator is passed:
when a scheduling backend must redraw a loser whose promised firing time
`tstar` overran, it draws from that loser's *own* keyed stream (see
[Randomness Ownership](randomness.md)), which is exactly what keeps a forced
firing coupled across two samplers built from the same seed. Second, there is
an **unbiasedness precondition**: the keep-if-later repair a scheduling
backend applies (a survivor whose stored schedule already exceeds `tstar`
keeps it) is proven correct when `tstar` is the current race's own decision
time, which is how a branching estimator always calls it. Choosing `tstar` by
peeking at a survivor's stored schedule biases the kept branch. See
[Contract and Invariants](contract.md).

## Capability traits

Not every sampler can answer every question. The exponential-only samplers
re-read a memoryless rate from the state at each step and keep no per-clock
enabling times, so age is undefined for them; only the backend that retains
draw randomness can hand back a retained uniform. Four Holy-style traits
(functions of the sampler *type*, so they resolve at compile time) let an
estimator check support at construction:

- [`supports_enabled_ages`](@ref) — the sampler keeps a per-clock
  enabling-time table and answers `enabled_ages`.
- [`supports_force`](@ref) — the sampler implements `force_fire!`.
- [`supports_retained_draw`](@ref) — the sampler can report the
  survival-space uniform behind a scheduled draw.
- [`supports_carry`](@ref) — the sampler can re-evaluate a still-enabled
  clock's distribution *deterministically*, consuming no fresh randomness;
  only such a sampler may be constructed with `coupling=:carry`
  (see [Re-evaluation and the Sampler's Coupling](reenable.md)).

The support matrix for the samplers in this package:

| sampler | `enabled_ages` | `force_fire!` | `retained_draw` | `carry` |
|:--------|:--------------:|:-------------:|:---------------:|:-------:|
| `CombinedNextReaction` | yes | yes | yes | yes |
| `FirstToFire` | yes | yes | no | yes |
| `FirstReaction` | yes | yes | no | no |
| `DirectCall`, `MultipleDirect`, `RSSA`, `PSSACR`, `Petri` | no | no | no | no |

`FirstToFire` can carry without a retained uniform because it reconstructs
the draw from its stored absolute firing time; `FirstReaction` redraws every
clock at every `next`, so it has no retained schedule to move but also
nothing to get wrong when forcing (forcing is just clock removal). The
exponential-only samplers draw for the race as a whole rather than for
individual clocks, so none of the aged verbs apply to them.

A verb called on a sampler that does not support it throws an
`ArgumentError` naming the missing trait — a clear failure at the estimator's
boundary rather than a `MethodError` deep inside a backend:

```julia
using CompetingClocks: DirectCall, enabled_ages
dc = DirectCall{Symbol,Float64}()
enabled_ages(dc, 0.0)
# ERROR: ArgumentError: enabled_ages is not supported by DirectCall{…};
# supports_enabled_ages(::DirectCall{…}) is false. An exponential-only
# sampler keeps no enabling-time table, so age is undefined for it.
```

An estimator should test the traits once, at its own construction, and refuse
a sampler that lacks what it needs.

## Which estimator needs which verb

The sampler design document (`sampler_design.tex` in the WorldTimer research
notes) measured, estimator by estimator, which verbs each one actually
consumed. The result is a strict tier structure:

- **Score function** (likelihood-ratio) estimators consume *nothing beyond a
  correct trajectory record*: the four simulation verbs plus the
  [`TrajectoryRecorder`](@ref)'s firing sequence. They re-weight recorded
  paths; they never touch a live sampler.
- **Pathwise / IPA** (infinitesimal perturbation analysis — differentiate the
  outcome along a fixed realization of the randomness) adds exactly one item:
  the retained draw `(te, u)` of each firing, which the record already
  carries. It replays firing times through the retained-draw identity at a
  perturbed parameter. For a *live* IPA run, a sampler constructed with the
  `:carry` re-evaluation coupling (the default for `CombinedNextReaction` and
  `FirstToFire`; see [`reenable!`](@ref)) is the sampler-side requirement.
- **Weak-derivative / branching** estimators add the cluster of live verbs:
  `enabled_ages` to build the selection law at the branch point, `clone` to
  copy the running sampler, `force_fire!` to impose the branched firing, and
  `rekey_streams!` to give a branch its own randomness afterward (see
  [Randomness Ownership](randomness.md)).

Every estimator, in addition, rebuilds distributions at the parameter of
interest on *its own* side of the boundary — the sampler itself only ever
sees plain `Float64` distributions (the "primal boundary" described in
[Contract and Invariants](contract.md)).

## Reference

```@docs
supports_enabled_ages
supports_force
supports_retained_draw
supports_carry
enabling_times
enabled_ages
retained_draw
force_fire!
```
