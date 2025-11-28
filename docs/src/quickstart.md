# Quickstart

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
sampler = SamplingContext(ClockKey, Float64, rng)
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
    - `fire!` that event, which also disables the event.
    - Handle changes to state and enabled events.

You tell it what *could* happen with `enable!()`, ask what happens `next()`,
and decide what you want to `fire!()`. Why is `next()` separate from `fire!()`?
So that you can stop a simulation at a fixed time rather than after an event.

## Next Steps

Learn More Techniques

 - State management: [Sample Main Loop](@ref "Sample Main Loop") - Better state handling.
 - Choosing samplers: [Choosing a Sampler](@ref "Choosing a Sampler") - When to use which algorithm.
 - Real examples:
   - [SIR Model](sir.md)
   - [Reliability](reliability.md)
   - [Birth-Death Process](constant_birth.md)
   - [Gene Expression](gene_expression.md)

Advanced Features

 - [Get likelihoods](@ref "Importance Sampling for Simulation")
 - [Variance reduction](@ref "Common Random Numbers")
 - [Build a framework](@ref "Integration Guide")

API Reference

 - [Build a Context](@ref "Building a Context")
 - [`CompetingClocks.enable!`](@ref) - Parameter documenation.
 - [`CompetingClocks.next`](@ref) - Return value details.
 - [`CompetingClocks.SamplerBuilder`](@ref) - Configuration options.
