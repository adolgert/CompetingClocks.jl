```@meta
CurrentModule = CompetingClocks
```

# Migration Guide: 0.3 to 0.4

Version 0.4 inverts the ownership of randomness: samplers own per-clock keyed
streams, and no sampler verb takes a random number generator. If your code
uses only the high-level [`SamplingContext`](@ref) — `enable!(ctx, ...)`,
`next(ctx)`, `fire!(ctx, ...)` — it is likely to run unchanged, because the
context API never threaded an rng through its verbs. Every break below is at
the low-level sampler layer or in the common-random-numbers machinery. The
concepts behind the changes are explained in
[Randomness Ownership](randomness.md).

## The `rng` argument is gone from sampler verbs

`enable!`, `next`, and `jitter!` on a low-level sampler no longer accept a
random number generator. The sampler draws from its own keyed streams.

Before (0.3):

```julia
rng = Xoshiro(2947223)
sampler = FirstToFire{Int,Float64}()
enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
when, which = next(sampler, 0.0, rng)
jitter!(sampler, when, rng)
```

After (0.4):

```julia
sampler = FirstToFire{Int,Float64}(2947223)   # the seed moves to the constructor
enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
when, which = next(sampler, 0.0)
jitter!(sampler, when)
```

If you kept an rng only to feed the sampler, delete it. If your *model* draws
its own random numbers (initial conditions, state-dependent choices), keep
your rng for those — the inversion applies to sampler verbs only.

## Sampler constructors take a seed

Every low-level sampler constructor accepts a seed that selects its family of
per-clock streams, with a fixed default so a bare constructor stays
deterministic:

```julia
CombinedNextReaction{Symbol,Float64}(0x5EED5)
FirstToFire{Int,Float64}(2947223)
DirectCall{Int,Float64}()               # default seed
```

Replications that used to come from advancing one shared rng now come from
distinct seeds — one sampler (or one re-key) per replication:

```julia
for rep in 1:n
    s = CombinedNextReaction{Symbol,Float64}(base_seed + rep)
    # ... run ...
end
```

A `SamplingContext` still takes an rng at construction; it uses it to draw
the sampler's stream seed (and to seed split copies), so passing equal-seeded
rngs to two contexts still controls reproducibility.

## `clone` changed meaning; `similar_sampler` is the old `clone`

In 0.3, `clone(sampler)` returned an *empty* sampler of the same type. In
0.4, [`clone`](@ref) is a **full-state copy** — clock state, schedules, and
the keyed streams with every generator's state and occurrence counts — so the
clone and the original run identical trajectories until one is re-keyed. The
old empty-copy behavior moved to [`similar_sampler`](@ref).

Before (0.3):

```julia
fresh = clone(sampler)              # empty sampler, same type and options
```

After (0.4):

```julia
fresh   = similar_sampler(sampler)  # empty sampler, same type and options
coupled = clone(sampler)            # full copy: runs the SAME future
rekey_streams!(coupled, newseed)    # ...until you re-key it
```

Audit every `clone(sampler)` call site: if it fed an array of independent
runs, it must become `similar_sampler` (or a `clone` followed by
`rekey_streams!`), otherwise your "independent" runs will be perfectly
coupled. The context-level `clone(ctx, rng)` keeps its old meaning — a fresh,
empty context built like the original — and re-keys the new sampler from the
rng you pass.

## `CommonRandom`, `freeze_crn!`, `reset_crn!`, `misses`/`misscount` are removed

The 0.3 common-random-numbers machinery recorded draws in a warm-up phase and
replayed them, with a miss count for draws that had no recording. All of it
is gone: the `CommonRandom` decorator, the `common_random=true` builder flag,
`freeze_crn!`, `reset_crn!`, and the miss counters. Coupling now falls out of
stream ownership, with no warm-up and no misses.

Before (0.3):

```julia
builder = SamplerBuilder(KeyType, Float64; common_random=true)
ctx = SamplingContext(builder, rng)
for i in 1:warm_up
    reset!(ctx)
    run_simulation(model, ctx)
end
freeze_crn!(ctx)                   # lock the recorded draws
for i in 1:draw_cnt
    reset!(ctx)
    run_simulation(perturbed_model, ctx)  # replays recorded draws
end
```

After (0.4):

```julia
# Same seed ⇒ coupled paths. No builder flag, no warm-up, no freeze.
base_ctx = SamplingContext(SamplerBuilder(KeyType, Float64), Xoshiro(trial))
pert_ctx = SamplingContext(SamplerBuilder(KeyType, Float64), Xoshiro(trial))
baseline  = run_simulation(model, base_ctx)
perturbed = run_simulation(perturbed_model, pert_ctx)
```

Each clock's `k`-th draw is the same variate in both runs wherever the two
trajectories enable the same clocks, which is strictly stronger than the old
call-order replay: inserting an extra event in one run no longer decouples
everything after it. See the [CRN recipe](randomness.md) for the low-level
form (`clone` + `rekey_streams!`) and
[Common Random Numbers](commonrandom.md) for a worked context-level example.

Passing `common_random=true` to `SamplerBuilder` is now an error (the keyword
no longer exists); remove it.

## `split!` no longer takes fresh rngs from you

Splitting a context into divergent copies used to rely on jittering with a
fresh rng per copy. Now [`split!`](@ref) copies clocks *and* stream states
faithfully, then re-keys each copy's streams from that copy's context rng and
jitters off the fresh streams. If you implemented splitting by hand, the
sequence is:

```julia
copy_clocks!(dst, src)                              # faithful: coupled
rekey_streams!(dst.sampler, rand(dst.rng, UInt64))  # decouple
jitter!(dst.sampler, time(src))                     # resample enabled clocks
```

Jittering *without* the re-key reproduces the original's draws exactly — a
faithful copy re-draws from identical generator states — so the re-key is not
optional.

## Re-enabling an enabled clock: prefer the explicit `reenable!`

This is a behavioral footnote rather than a break. Calling `enable!` on an
already-enabled clock still re-evaluates its distribution in place, but which
pathwise coupling that implements is chosen silently by the backend:
`CombinedNextReaction` carries its retained draw through the change, while
`FirstToFire`, `FirstReaction`, and the exponential-only samplers redraw. The
couplings agree in law, so no statistical result changes — but if you are
coupling runs or differentiating along paths, the silent difference matters.
0.4 adds [`reenable!`](@ref)`(..., :carry | :redraw)` so the choice is yours
per call; new code should use it, and code that relied on
`CombinedNextReaction`'s silent carry should say `:carry` explicitly rather
than depend on the backend default. See
[Re-evaluation Couplings](reenable.md).

## New in 0.4, not breaking

Alongside the breaking changes, 0.4 adds the estimator-facing surface — the
[`TrajectoryRecorder`](@ref) ([Recording Trajectories](recording.md)), the
[`enabled_ages`](@ref)/[`retained_draw`](@ref)/[`force_fire!`](@ref) verbs
with their capability traits ([Estimator-Facing Verbs](estimator_verbs.md)),
and the invariants they rest on ([Contract and Invariants](contract.md)).
None of these change existing simulation code.
