# Integration Guide

How your code can use CompetingClocks.jl as a sampling engine.

## Core Integration Pattern

CompetingClocks manages _when_ events happen. Your framework manages _what_ happens.

Your framework provides:
 * State representation
 * Event-handling logic
 * A user-facing API

CompetingClocks provides:
 * Event timing and scheduling
 * Distribution sampling
 * Variance reduction and likelihood calculation

A minimal integration loop:
```julia
function run_simulation(model, sampler, rng, end_time)
    initialize_events!(model, sampler)  # Your code

    when, which = next(sampler)
    while !isnothing(which) && when < end_time
        fire!(sampler, which, when)
        handle_event!(model, sampler, which, when, rng)  # Your code
        when, which = next(sampler)
    end
end
```

## Design Decisions

### Clock Key Type

How will you identify events? CompetingClocks.jl indexes events using an
immutable clock key that you provide. Examples:

 * `(:infect, infectious, susceptible)` - `Tuple{Symbol,Int,Int}` or `Tuple`
 * `392` - `Int64`
 * `:i => 37` - A Pair for use with `Gen.jl`.

A custom struct type for events works well. Using a concrete type is more
performant for sampling.

### Threading

This package doesn't use Tasks or thread guards because it is so common to
run multiple copies of a simulation, each on its own thread. If you run
multiple simulations, it can be more correct to use random number generators
that work in parallel, such as those in the `Random123` package or the `jump!()`
method for the `Xoshiro` sampler in the `Random` package.

### State and Dependencies

Your simulation handles all state. `CompetingClocks` only knows about enabled
events. It can be very helpful for performance if you use information about
relationships among events to limit updates to the samplers. There are a few
kinds of dependency graphs that may help.

 * Dependency Graph---When one event fires, it changes state. Which other events could be
   enabled or disabled because of the changed state?
 * Rate Dependency Graph---When one event fires, it changes state. Which other events have transition
   *rates* that depend on that state? The event may remain enabled, but its
   distribution of firing times may have a different parameter value.
 * Event-to-event Graph---When one event fires, you might know directly which other events have
   enabling rules or transition rates that depend on that first event, without
   needing to ask what state changed and what events depend on the state.

How you represent a dependency graph depends on the problem.

## Common Integration Patterns

### Agent-based Models

The clock key will identify an action and an agent, along with either
another agent or a location.

```julia
struct ClockKey
    action::Symbol
    who::Int
    other::Int # other agent or location
end
```

It's also easy to use a `Tuple` as a ClockKey until you figure out the space
of all possible keys.

## Spatial Simulation

Split the sampler into groups so that each region is faster.

```julia
builder = SamplerBuilder(EventKey, Float64)
for region in spatial_grid
    add_group!(builder, region.name => (key, dist) -> key.location == region.id)
end
```

### Chemical Reaction Networks
Map reactions to events.
```julia
function enable_reaction!(sampler, reaction, state, params)
    propensity = calculate_propensity(reaction, state, params)
    enable!(sampler, reaction.id, Exponential(1/propensity))
end
```

## Advanced Features

There is extra work in a simulation to support some of the sampler features.

### Exposing likelihood calculation

The step log-likelihood and path log-likelihood can be found from the sampler
at any point during a simulation. They are always calulated relative to the
last event to `fire!()`, even if `next()` has been called.
```julia
builder = SamplerBuilder(KeyType, Float64;
    path_likelihood=true)
sampler = SamplingContext(builder, rng)
# After simulation, the log_prob is a Float64.
log_prob = pathloglikelihood(sampler, end_time)
```
The resulting value is the likelihood of a single trajectory of the system.
The `pathloglikelihood` is a regular function and can be differentiated.

## Variance reduction

Common random numbers (CRN) need no builder flag and no warm-up phase: every
sampler owns per-clock keyed random streams, and two contexts built from the
same seed consume identical randomness per clock. To compare a baseline model
against a perturbed one, build both contexts from equal-seeded rngs:

```julia
for trial in 1:trial_cnt
    base_ctx = SamplingContext(SamplerBuilder(KeyType, Float64), Xoshiro(trial))
    pert_ctx = SamplingContext(SamplerBuilder(KeyType, Float64), Xoshiro(trial))
    baseline[trial]  = run_simulation(model, base_ctx)
    perturbed[trial] = run_simulation(perturbed_model, pert_ctx)
end
```

The paired difference `perturbed .- baseline` has far smaller variance than
independent seeds would give, because each clock's draws are shared between
the pair even when the two runs fire events in different orders. Each
simulation copy still needs its own copy of the model state. For the
lower-level primitives — a full-state `clone` for a coupled continuation and
`rekey_streams!` to decouple it — see [Randomness Ownership](randomness.md)
and the [Common Random Numbers](commonrandom.md) example. If you migrated
from 0.3's `common_random=true` recorder, the [Migration Guide](migration_04.md)
maps each removed function to its replacement.

### Path splitting for rare events

This is a simple form of importance sampling. When a simulation gets near
a space of interesting results, take the current trajectory and make multiple
copies of the simulation, both the state and the samplers. Weight the results
of each copy by one over the number of copies.

CompetingClocks offers a `split!` function that copies the state to a vector
of clocks and keeps a `split_weight` as a record of how many times splitting
was done. Take a look at the implementation of [`CompetingClocks.split!()`](@ref).

Your code might keep the sampler in the main simulation framework, in which case
you would include the sampler in the cloning and copying process.

## Error Handling and Debugging

If you call the sampler builder with `debug=true` then there will be logging.
```julia
builder = SamplerBuilder(KeyType, Float64; debug=true)
sampler = SamplingContext(builder, rng)
```
Then run with logging on.
```julia
with_logger(ConsoleLogger(stdout, Logging.Debug)) do
    observed, importance = run_epochs(100)
end
```

There is also a [chaos-monkey](https://en.wikipedia.org/wiki/Chaos_engineering) option called a `Petri` sampler.
This sampler always advances time by the same time step `dt=1.0` by default,
and it ignores the clock distributions. It chooses any enabled sampler at
each step, with the goal of finding unusual sampling paths. This is surprisingly
effective at finding logic errors in code to change state.

## Performance Considerations

**Type Stability**---Use a concrete type for the
ClockKey for more speed.

**Sampler Choice**---If there are a lot of clocks enabled at the same time,
the the sampler choice can be important. Experts in traffic simulation use
hierarchical samplers to split sampling by geographic location but also
split fast events from slow events. There are also fine-grained choices for
samplers about whether clock keys are reused.

**Memory Allocation**---Instead of making a new sampler for each run, use
`reset!(sampler)` between runs to clear out the clocks.
