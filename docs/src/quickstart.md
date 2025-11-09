# CompetingClocks.jl Quickstart

Learn CompetingClocks.jl basics by building a minimal simulation.

## The Simulation

A birth-death process, which is a population that grows and shrinks stochastically.
 * Birth events at rate 2.0 per time unit per individual.
 * Death events are Gamma-distributed for each individual.

```julia
using CompetingClocks
using Distributions
using Random

rng = Xoshiro(9834223)
person_cnt = 10
population = Set(1:person_cnt)

# Events will be identified by a tuple of type and individual.
const ClockKey = Tuple{Symbol,Int}
# Create a sampler.
builder = SamplerBuilder(ClockKey, Float64)
sampler = SamplingContext(builder, rng)
enable!(sampler, (:birth, 0), Exponential(inv(length(population) * 2.0)))
for mort in population
    enable!(sampler, (:death, mort), Gamma(inv(1.5)))
end
when, event = next(sampler)
while !isnothing(event) && when < 100.0
    fire!(sampler, event, when)
    if event[1] == :birth
        person_cnt += 1
        push!(population, person_cnt)
        enable!(sampler, (:death, person_cnt), Gamma(inv(1.5)))
        enable!(sampler, (:birth, 0), Exponential(inv(length(population) * 2.0)))
    elseif event[1] == :death
        delete!(population, event[2])
        enable!(sampler, (:birth, 0), Exponential(inv(length(population) * 2.0)))
    else
        error("Unknown event $event")
    end
    when, event = next(sampler)
end
```

## The Core Pattern

 1. Create a state.
 2. Define a type for events.
 3. Create a sampler.
 4. Enable initial events.
 5. Ask what happens next and
    - `fire!` that event.
    - Handle changes to state and enabled events.

You tell it what *could* happen with `enable!()`, ask what happens `next()`,
and decide what you want to `fire!()`. Why is `next()` separate from `fire!()`?
So that you can stop a simulation at a fixed time rather than after an event.

## Next Steps

Learn More Techniques

 - State management: your-first-sim.md - Better state handling.
 - Choosing samplers: choosign-samplers.md - When to use which algorithm.
 - Real examples:
   - SIR
   - reliability

Advanced Features

 - Get likelihoods: likelihood-tutorial
 - Variance reduction: variance-reduction.md
 - Build a framework: integration-guide.md

API Reference

 - api.md#enable - Full parameter documentation
 - api.md#next - Return value details
 - api.md#samplerbuilder - Configuration options
