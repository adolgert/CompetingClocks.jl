```@meta
CurrentModule = CompetingClocks
```

# Recording Trajectories

A derivative estimator, a likelihood replay, or an offline analysis needs the
same raw material: the ordered sequence of clock firings that a simulation
produced, together with enough information to reconstruct the probability of
each firing. This page shows how to capture that record with a
`TrajectoryRecorder`, what each entry of the record contains, and the identity
that every recorded firing satisfies. The estimators that consume the record
are described in [Estimator-Facing Verbs](estimator_verbs.md).

## Attaching a recorder

A `TrajectoryRecorder` is a *watcher*: an observer that receives the same
`enable!`, `disable!`, and `fire!` notifications that the sampler receives,
and whose product is the ordered list of firings. Attach one to a running
[`SamplingContext`](@ref) with [`with_recorder`](@ref), which returns an
extended context and the recorder itself. Drive the extended context from then
on; read the firings from the recorder.

```julia
using CompetingClocks
using CompetingClocks: with_recorder, close_record!, recorded_firings
using Distributions
using Random: Xoshiro

ctx = SamplingContext(SamplerBuilder(Symbol, Float64; method=NextReactionMethod()),
                      Xoshiro(20260709))
ctx, rec = with_recorder(ctx)

# ... run the simulation on ctx ...

close_record!(rec, 10.0)       # observation ran to horizon 10.0
recorded_firings(rec)          # Vector{ClockFiredRecord}
```

Call [`close_record!`](@ref) when the run ends. Its `horizon` argument stamps
the time the trajectory was observed to, which is where the survival
(censoring) terms of every clock still enabled at the end of the run extend
to. A likelihood replay needs the horizon to account for the clocks that were
racing but never fired.

## What the record contains

Each firing is stored as a [`ClockFiredRecord`](@ref) with six fields:

| field | meaning |
|:------|:--------|
| `clock` | the key of the clock that fired |
| `when` | the absolute firing time |
| `te` | the enabling time — the zero-point of the clock's lifetime distribution, in absolute time |
| `distribution` | the lifetime distribution that was current at the firing |
| `u` | the survival-space uniform of the clock's total lifetime, `ccdf(distribution, when - te)` |
| `logu` | `log(u)`, the log-survival coordinate |

The first four fields say *what happened*. The last two say *which random
draw produced it*, and they are the part that makes the record useful to
derivative estimators.

## The retained-draw identity

Every recorded firing satisfies

```
when == te + invlogccdf(distribution, logu)      equivalently      u == ccdf(distribution, when - te)
```

In plain English: `u` is the probability that this clock, with this lifetime
distribution measured from its enabling time `te`, would have survived past
the moment it actually fired. A uniform random number and a firing time are
two coordinates for the same draw — given the distribution and `te`, either
one determines the other. `u` lives in *survival space* (a survival
probability near 1 means an early firing; near 0 means a late one), so it is
the natural coordinate for coupling and for pathwise derivatives: perturb the
distribution's parameters while holding `u` fixed, and the identity tells you
exactly how the firing time moves.

The identity is stated — and stored — in log space. `logu = log(u)` is kept
alongside `u` because a firing deep in a distribution's tail underflows `u`
to `0.0` in 64-bit floating point while `logu` stays finite, and it is `logu`
(not `u`) that inverts back to `when` without loss of precision.

Two details are worth knowing before you consume the record.

First, the recorder is **sampler-agnostic**. It never asks the sampler for
its internal random draw. It knows `(distribution, te)` from the `enable!`
notification it received, and at fire time it back-calculates
`logu = logccdf(distribution, when - te)` from the observed firing time. A
first-reaction sampler that never draws by inversion and a next-reaction
sampler that stores a survival uniform internally therefore produce identical
records for the same trajectory.

Second, `u` is the uniform of the clock's **total lifetime**, not a uniform
conditioned on survival to the enable call. When a clock is enabled with its
enabling time shifted into the past (`te < when` at the enable, the
left-shift idiom of [Re-enabling and Memory](memory_idioms.md)), the identity
stays anchored at `te`, so a clock that was already aged when the sampler
first saw it still records the survival of its whole life. The two uniforms
differ — the conditioned one would be
`ccdf(dist, when - te) / ccdf(dist, age)` — and replaying the wrong one gives
the wrong likelihood.

## A worked example: machines that break and get repaired

Five machines each carry a failure clock while working and a repair clock
while broken; firing one enables the other. This is the machine-repair
pattern from the package's own recorder test
(`test/test_recorder.jl`).

```julia
using CompetingClocks
using CompetingClocks: with_recorder, close_record!, recorded_firings
using Distributions
using Random: Xoshiro

function run_machine_repair!(ctx, nmachines, horizon, fail_dist, repair_dist)
    for i in 1:nmachines
        enable!(ctx, (:fail, i), fail_dist(i))
    end
    while true
        when, which = next(ctx)
        when > horizon && break
        fire!(ctx, which, when)
        kind, i = which
        if kind === :fail
            enable!(ctx, (:repair, i), repair_dist(i))
        else
            enable!(ctx, (:fail, i), fail_dist(i))
        end
    end
end

fail_dist(i) = Weibull(1.7, 1.0 + 0.1 * i)   # non-exponential, so age matters
repair_dist(i) = Exponential(0.5)
horizon = 40.0

ctx = SamplingContext(SamplerBuilder(Tuple{Symbol,Int}, Float64;
                                     method=NextReactionMethod()),
                      Xoshiro(20260709))
ctx, rec = with_recorder(ctx)
run_machine_repair!(ctx, 5, horizon, fail_dist, repair_dist)
close_record!(rec, horizon)

for fr in recorded_firings(rec)
    # Every firing satisfies the retained-draw identity.
    @assert isapprox(fr.when, fr.te + invlogccdf(fr.distribution, fr.logu); atol=1e-9)
end
```

Because the identity is a contract obligation on every firing (see
[Contract and Invariants](contract.md)), the loop's assertion holds on every
sampler in the package, and the recorded firing sequence is a *sufficient
statistic* for the path: a score-function estimator or an offline pathwise
replay can reconstruct the trajectory's likelihood — and its derivative in
the model's parameters — from the record alone, without touching the sampler
again.

## Recorder life cycle

- `reset!(recorder)` clears the firing history, the live clock table, and the
  horizon, so one recorder can serve many runs.
- [`isclosed`](@ref) and [`horizon`](@ref) report whether and where the
  record was closed; check `isclosed` to distinguish an unclosed record from
  one closed at time zero.
- `clone(recorder)` returns a fresh, empty recorder of the same type (a clone
  seeds a new run, it does not duplicate a history), while `copy_clocks!`
  carries both the live table and the firing history to a destination
  recorder, so a split or branched continuation inherits the path recorded so
  far.
- [`attach_watcher`](@ref) is the general form of `with_recorder`: it rebuilds
  the context around an extended watcher tuple, sharing the sampler and its
  state with the original. Use only the returned context afterward.

Forced firings — the [`force_fire!`](@ref) verb described in
[Estimator-Facing Verbs](estimator_verbs.md) — reach the recorder through the
same notification as ordinary firings, so a record that contains them still
satisfies the identity on every entry.

## Reference

```@docs
TrajectoryRecorder
ClockFiredRecord
with_recorder
attach_watcher
close_record!
recorded_firings
isclosed
horizon
```
