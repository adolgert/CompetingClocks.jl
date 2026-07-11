```@meta
CurrentModule = CompetingClocks
```

# Contract and Invariants

This page collects the obligations that hold across the sampler layer — the
guarantees an estimator may rely on and the preconditions a caller must
honor. The manual pages ([Recording Trajectories](recording.md),
[Estimator-Facing Verbs](estimator_verbs.md),
[Randomness Ownership](randomness.md),
[Re-evaluation and the Sampler's Coupling](reenable.md)) show how to use the machinery; this
page states what must always be true.

## The retained-draw identity

Every firing record produced by a [`TrajectoryRecorder`](@ref) satisfies

```
when == te + invlogccdf(distribution, logu)      with      u == exp(logu) == ccdf(distribution, when - te)
```

where `distribution` and `te` are the lifetime distribution and enabling time
current at the firing. This is a **contract obligation on every firing**, not
a property of one backend: the recorder back-calculates `u` from the observed
firing time, so any correct sampler produces records satisfying it, and any
record violating it indicates a broken enable/disable/fire sequence. The
identity is what makes the record a sufficient statistic for likelihood
replay and pathwise derivatives. It is stated in log-survival space because
`u` underflows in deep tails while `logu` stays finite and inverts without
loss.

The `u` in the identity is the uniform of the clock's **total lifetime**
measured from `te` — not a uniform conditioned on survival to the moment the
sampler was told about the clock. For a left-shifted enabling (`te < when` at
the enable call) the two differ, and the total-lifetime one is the recorded
coordinate.

## `fire!` versus `disable!`

The two verbs that end a clock's enabled episode act differently on the
clock's draw:

- `fire!` **realizes** the draw. The firing time is observed exactly; no
  residual randomness remains; a re-enabled clock must draw fresh.
- `disable!` **censors** the draw. The simulation learned only that the clock
  had not yet fired, so a sampler is entitled to retain the draw's remaining
  randomness and reuse it if the clock is enabled again.

The distinction is invisible for exponential clocks and load-bearing for
everything else. In this package, `CombinedNextReaction` is the sampler that
retains residual draw randomness across enabling episodes (the
Anderson/Gibson–Bruck reuse property); for every other sampler the two verbs
act identically on stored state. Callers must use the verb that matches what
happened: `disable!` on a clock that actually fired lets a reusing sampler
resurrect a consumed draw, and `fire!` on a clock that was merely switched
off discards randomness the sampler was entitled to keep — the law survives,
but draw reuse and common random numbers do not. `force_fire!` consumes the
chosen clock's draw as `fire!` does, not as `disable!` does.

## `next` is a reservation, and `when` is monotone

The `(when, which)` returned by `next` is a *reservation*, not a commitment.
It is valid only until the next `enable!`, `disable!`, or `fire!` changes the
sampler's state, and repeated calls without a state change need not return
identical values across samplers (`FirstReaction` redraws per call;
`CombinedNextReaction` returns its cached reservation). The two supported
responses are to fire it or to decline it and stop; there is deliberately no
`peek`.

The `when` argument passed to `next` must not decrease from one call to the
next, and must never advance past a pending firing time without that event
being fired. The time of the most recently fired event is always a safe
choice. Within this invariant, re-querying at a later `when` is legal — it is
how hierarchical samplers re-query losing sub-samplers at the advanced global
time — and all samplers agree in law. Outside it, they legitimately disagree.

## Deterministic enabled-key order

[`enabled_ages`](@ref) returns clocks **sorted by key**, and the order is
load-bearing, not cosmetic: a forced-selection probability-mass function (the
branch-selection step of a weak-derivative estimator) indexes into this order,
so two samplers, or a sampler and its clone, must present the same clocks at
the same positions. Any future sampler or model integration must preserve a
deterministic enabled-key order.

## The `force_fire!` unbiasedness precondition

`force_fire!(sampler, key, tstar)` repairs each surviving clock by
keep-if-later / redraw-if-passed: a stored schedule strictly beyond `tstar`
is kept exactly, and a schedule at or before `tstar` (a promise the force
overran) is redrawn from the lifetime conditioned on survival past `tstar`.
The keep branch is proven unbiased **when `tstar` is the current race's own
decision time** — the time at which the race would have been resolved — which
is how a branching estimator always calls it. Choosing `tstar` by peeking at
a survivor's stored schedule conditions on the very draw being kept and
biases the kept branch. Callers own this precondition; the sampler cannot
check it.

## Storable is not replayable

The number a sampler *stores* for a clock is generally not the coordinate an
estimator can *replay*. Two instances of this rule matter here:

1. The recorder stores the **total-lifetime** survival uniform, not the
   uniform conditioned on survival to the enable call (see
   [the identity](#The-retained-draw-identity) above). Replaying the
   conditional uniform against the unconditioned distribution gives the wrong
   time.
2. `CombinedNextReaction` stores each clock's survival in whichever sampling
   space suits the distribution family (linear or log), and — for a
   left-shifted enabling with `te < t0`, where `t0` is the last (re)schedule
   time — within the distribution *truncated* at the age `t0 - te`.
   [`retained_draw`](@ref) therefore cannot just read the stored number; it
   normalizes to the total-lifetime log-survival coordinate by adding the
   truncation shift:

   ```
   logu = log_survival(space, survival) + (te < t0 ? logccdf(dist, t0 - te) : 0)
   ```

   After that normalization, one identity —
   `tentative == te + invlogccdf(dist, logu)` — holds for every distribution
   family and every enabling shift.

The general statement: the **canonical coupling coordinate is the
survival-space uniform `u` (equivalently `logu`) of the total lifetime from
`te`**. Raw uniforms consumed from a generator are an implementation detail —
a single `rand(gen, dist)` may consume several native words — and are never
the replay or coupling coordinate.

## Stream addressing

The [`KeyedStreams`](@ref) family obeys three invariants:

- **Per-key generators are seeded once**, from `hash((seed, key))`, at the
  key's first draw. The hash selects a seed; it never synthesizes uniforms.
  A key's stream advances only when that key draws, so a key's draw sequence
  is independent of every other key's activity.
- **Coupling is by `(clock, occurrence)`** in the survival-space coordinate,
  as above — never by matching positions in a global consumption order.
- **Occurrence counts are clone state.** `copy(streams)` — and therefore
  [`clone`](@ref) of a sampler — carries generator *states and counts*, so a
  clone continues the original's per-key sequences exactly. Forgetting either
  would silently couple or decouple a copy's future; that omission is the bug
  class the design exists to prevent. [`rekey_streams!`](@ref) is the one
  sanctioned way to decouple, and [`split!`](@ref) composes copy + re-key +
  jitter.

The reserved race generator lives *inside* `KeyedStreams`, so copying and
re-keying handle it in the same place as the per-key generators.

A fourth requirement falls on the *caller*: **clock keys must hash by
content.** Because seeding uses `hash((seed, key))`, any key component whose
`Base.hash` falls back to `objectid` — notably `@enum` values — makes the
stream seed process-dependent, silently breaking same-seed reproducibility
and cross-module key identity. Symbols, numbers, strings, and tuples of them
are safe. If a key contains an enum `E`, define
`Base.hash(x::E, h::UInt) = hash(Symbol(x), h)` beside the enum. (This trap
was found in practice by the elevator example's direction enum.)

## The primal boundary (dual firewall)

Samplers run on plain `Float64` — always. When a caller's distribution
carries derivative-tracking number types (such as `ForwardDiff.Dual`), the
context strips them at the boundary via `primal_distribution`, handing the
sampler a value-only shadow while the differentiable copy stays on the
likelihood watcher. This is sometimes called the **dual firewall**: dual
numbers (value-plus-derivative pairs) never cross into sampler internals.
Derivatives enter only in estimator-layer computations — likelihood
accumulation and offline replays of the recorded trajectory — where the
estimator rebuilds distributions at its parameter of interest. The payoff is
that every sampler stays type-stable and compiles without any automatic
differentiation package, and no derivative information can silently alter
which trajectory is sampled.

## Capability traits fail fast

A sampler declares what it supports through the traits
[`supports_enabled_ages`](@ref), [`supports_force`](@ref),
[`supports_retained_draw`](@ref), and [`supports_carry`](@ref); the default
for a new sampler type is `false` for all of them (opt-out until deliberate
opt-in). Calling a verb on a sampler that opted out throws an
`ArgumentError` naming the missing trait at the verb boundary — never a
`MethodError` from inside a backend. An estimator should check the traits it
needs **at its own construction**, so an unsupported pairing fails before any
simulation work is done.
