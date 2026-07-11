```@meta
CurrentModule = CompetingClocks
```

# Randomness Ownership

As of version 0.4, no sampler verb takes a random number generator. Each
sampler *owns* its randomness as a family of per-clock keyed streams, and a
random draw belongs to a clock — not to a position in a global call sequence.
This page explains the mechanism, the two kinds of sampler copies it makes
possible, and the recipe for common random numbers (CRN — the variance
reduction technique of running two parameter settings on the same underlying
randomness). If you are migrating code that passed an `rng` to `enable!` or
`next`, see the [Migration Guide](migration_04.md).

## Draws belong to clocks

A [`KeyedStreams`](@ref) holds one pseudo-random generator per clock key,
each seeded once from the pair `(seed, key)`, plus each key's *occurrence
count* — how many times that key has drawn. When clock `:a` needs a firing
time, the sampler draws from `:a`'s own generator; the draw is addressed by
the coordinate `(clock, occurrence)`.

The consequence is order independence: two samplers built from the **same
seed** consume identical randomness per clock *even when their events fire in
a different order*, because each clock's generator advances only when that
clock draws. Under the old convention — one global generator threaded through
every call — inserting a single extra event anywhere in a run shifted every
subsequent draw and silently decoupled everything after it.

```julia
using CompetingClocks: FirstToFire, enable!, next, fire!
using Distributions

dists = Dict(:a => Weibull(1.5, 1.0), :b => Exponential(1.3),
             :c => Gamma(2.0, 0.7), :d => Weibull(2.0, 0.9))

s1 = FirstToFire{Symbol,Float64}(0xC0FFEE)
for k in (:a, :b, :c, :d)
    enable!(s1, k, dists[k], 0.0, 0.0)
end

s2 = FirstToFire{Symbol,Float64}(0xC0FFEE)
for k in (:d, :c, :b, :a)          # the reverse enable order
    enable!(s2, k, dists[k], 0.0, 0.0)
end

@assert all(s1[k] == s2[k] for k in (:a, :b, :c, :d))  # identical schedules
```

Every sampler constructor accepts a seed —
`FirstToFire{Symbol,Float64}(0xC0FFEE)`,
`CombinedNextReaction{Int,Float64}(42)`, and so on — which selects the whole
family of streams. A sampler constructed without a seed uses a fixed default,
so a bare constructor is still deterministic. A [`SamplingContext`](@ref)
draws its sampler's stream seed once from the `rng` you pass at construction;
that is now the context rng's main job, so two contexts built from
equal-seeded rngs are coupled the same way two equal-seeded samplers are.

Some samplers (`DirectCall`, `MultipleDirect`, `RSSA`, `PSSACR`, `Petri`)
draw for the *race* as a whole — the time to the next event and which clock
wins — rather than for one clock at a time. Those draws come from a reserved
race stream ([`race_stream`](@ref)) carried inside the same `KeyedStreams`,
so copying and re-keying handle it together with the per-key generators. For
these samplers the coupling unit is the whole race, not the individual clock,
which is why per-clock CRN is most effective on the aged backends.

## Why no verb takes an RNG

Removing the RNG parameter is not a convenience; it is what makes coupling
*compositional*. When randomness is an argument, any caller can pass a
different generator and silently diverge two runs that were meant to be
coupled — or forget to, and silently couple two runs that were meant to be
independent. When randomness is owned and keyed, divergence is an explicit
act with a name (`rekey_streams!`), and the type system has no channel
through which call order can leak into the draws. The estimators that branch
and clone simulations (see [Estimator-Facing Verbs](estimator_verbs.md))
depend on exactly this property: a forced firing that must redraw a loser
draws from that loser's own stream, so a coupled clone pair redraws
identically.

## Coupled copies and fresh copies

Two verbs copy a sampler, with opposite intents:

- [`clone(sampler)`](@ref clone) is the **full-state** copy: clock state,
  distributions, schedules, *and* the keyed streams including every live
  generator's state and occurrence count. Running the clone and the original
  forward produces identical firing sequences — they are coupled. This is the
  primitive behind clone coupling and CRN.
- [`similar_sampler(sampler)`](@ref similar_sampler) is the **empty** copy:
  same type, same constructor options (including the seed), but no clock
  data and fresh streams. Use it to initialize an array of samplers or to
  start a new run.

To *decouple* a clone — give it a future of its own — re-seed its streams:

```julia
c = clone(s)                 # coupled: same draws ahead
rekey_streams!(c, 0x9E3779B9) # decoupled: fresh, independent stream family
```

[`rekey_streams!`](@ref) forgets every live generator and count and re-seeds
the family, so "pass a different RNG to diverge" is replaced by "re-key to
diverge". Note the interaction with `jitter!` (resampling the enabled clocks
at a stopping time): jittering a *faithful* clone reproduces the original's
draws, because the jitter re-draws from identical generator states. Divergence
requires the re-key first.

## `split!`: many divergent continuations

[`split!(dst, src)`](@ref split!) is the packaged form of that recipe for
path splitting (branching a promising run into several weighted
continuations). For each destination context it performs three steps: a
faithful `copy_clocks!` (clocks *and* stream states), a `rekey_streams!` from
a fresh seed, and a `jitter!` to resample the currently-enabled clocks off
the freshly-seeded streams. Each copy's `split_weight` is divided by the
number of copies so the branches stay properly weighted.

## The common-random-numbers recipe

To compare a baseline parameter `θ` against a perturbed `θ + h` with common
random numbers, build the two runs from the same seed. There is no recorder
to warm up, no replay mode, and no miss count — the coupling falls out of
stream ownership. The following pairs a machine-repair model at two parameter
values (this is the experiment from `test/test_crn_streams.jl`, where the
paired difference showed roughly fifty-fold variance reduction over
independent seeds):

```julia
using CompetingClocks: CombinedNextReaction, enable!, next, fire!
using Distributions

# Count repair completions by horizon T for break-time scale `theta`.
function repairs_by_horizon(seed, theta)
    s = CombinedNextReaction{Tuple{Int,Symbol},Float64}(seed)
    for i in 1:5
        enable!(s, (i, :break), Exponential(theta), 0.0, 0.0)
    end
    t, repairs = 0.0, 0
    while true
        when, which = next(s, t)
        (which === nothing || when > 50.0) && break
        i, kind = which
        fire!(s, which, when)
        if kind === :break
            enable!(s, (i, :repair), Exponential(0.5), when, when)
        else
            repairs += 1
            enable!(s, (i, :break), Exponential(theta), when, when)
        end
        t = when
    end
    return repairs
end

# CRN: the SAME seed drives both parameter settings of each pair.
diffs = [repairs_by_horizon(r, 1.05) - repairs_by_horizon(r, 1.0) for r in 1:1000]
```

Because machine `i`'s `k`-th break draw is the same underlying variate in
both runs of a pair, `diffs` estimates the parameter effect with far less
variance than independent seeds would — which means far fewer replications
for the same confidence. The same recipe works at the context level: build
two `SamplingContext`s from equal-seeded rngs, as shown in
[Common Random Numbers](commonrandom.md).

One caveat inherited from all CRN methods: the coupling helps only while the
two trajectories take similar events. If the parameter change redirects the
path enough that different clocks race, fewer draws are shared and the
variance reduction fades. Measure it — compare the paired variance against
independent seeds — before relying on it.

## Reference

The stream type and its verbs are documented in the
[API reference](reference.md); the two below are the low-level accessors a
custom sampler implementation would use.

```@docs
stream_for!
race_stream
```
