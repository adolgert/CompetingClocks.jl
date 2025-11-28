# Choosing a Sampler

How to make a sampler for events in CompetingClocks.jl.

## High-level Interface

Every event in a simulation has an identifying key which can be any immutable Julia
type. A `Tuple` key type works fine, but it will be more performant to use a
concrete type for the `KeyType`. The second argument
is a type to use for time. Basic usage:
```julia
sampler = SamplingContext(KeyType, Float64, rng)
```
The sampler takes a `Random.AbstractRNG` as its second argument.

Specify features for the sampler with keyword arguments.
These features determine the type of the sampler.

 * `method=nothing`---Setting this chooses a single sampling algorithm for all clocks.
 * `step_likelihood=false`---Setting this to `true` let's you calculate the
   the likelihood of the next event and time before firing it.
 * `path_likelihood=false`---Set to true in order to calculate likelihood of a
   whole trajectory at any point in the simulation. This is more efficient than
   adding up steps in the likelihood along the way.
 * `likelihood_cnt=1`---Applies when likelihoods are enabled and supports
   importance sampling. Specifies how
   many event distributions will be used to calculate a vector of likelihoods.
 * `common_random=false`---For variance reduction, turn on recording of
   random number usage during sampling.
 * `start_time=0`---Sometimes you want a simulation to start at a different time.
 * `debug=false`---Whether to print debug messages using the `Logging` package.
 * `recording=false`---This will create a vector of every enable and disable
   event for test and debug.

If you want to try a particular sampler instead of having the `SamplingContext`
pick one for you, use the `method` parameter.

```julia
sampler = SamplingContext(KeyType, Float64, rng; method=DirectMethod(:keep, :tree))
```

## Quick Reference

These are all of the samplers.

| Name           | Distribution support | Best use cases                 |
|----------------|----------------------|--------------------------------|
| [`FirstToFireMethod`](@ref)  | All                  | Fastest for simple simulation  |
| [`NextReactionMethod`](@ref)  | All                  | Best for common random numbers |
| [`FirstReactionMethod`](@ref) | All                  | Fastest for very few events    |
| [`DirectMethod`](@ref)         | Exponential-only     | Quicker likelihood calculation |
| [`RejectionMethod`](@ref)         | Exponential-only     | RSSA for large systems |
| [`PartialPropensityMethod`](@ref)         | Exponential-only     | For reaction networks. |
| [`PetriMethod`](@ref)          | All                  | Debug rare errors              |


## Hierarchical Samplers

For spatial simulations, chemical simulations that aren't well-mixed, or
simulations with some fast and many slow clocks,
it can help to split a sampler into multiple buckets so that each bucket can
update its own list of what fires next. This package does this by adding
sampler groups.

A sampler group is

 * A Symbol name for the sampler.
 * An inclusion function from `(ClockKey, Distribution)` to `Bool` that decides
   whether a given event belongs to this sampler.
 * An optional `method` to say what kind of sampler this group should use.

This example sends events with keys like `(:infect, 3)` tto the sampler called
`:forthright` and sends keys like `(:recover, 7)` to the sampler called
`:sparky`.
```julia
const KeyType = Tuple{Symbol,Integer}
builder = SamplerBuilder(KeyType, Float64)
add_group!(builder, :sparky => (x,d) -> x[1] == :recover, method=NextReaction())
add_group!(builder, :forthright=>(x,d) -> x[1] == :infect)
sampler = SamplingContext(builder, rng)
```

The resulting sampler will be hierarchical. For every event, it will choose the
soonest time among the soonest times from all sampler groups.
