# Delayed Clocks

## Introduction

There are common uses for delayed clocks, also called delayed reactions.

  1. State-independent uninterruptible delay: a susceptible-infectious-recovered (SIR) model where the I to R transition occurs after some fixed duration $\tau$.
  2. State-independent delay with state-independent interruption: consider the SIR model of type 1, but where there is an additional transition from I to death, with constant per-capita hazard $\mu$. Then the I to R transition may be interrupted by death with probability $1-e^{-\mu\tau}$, which is known at initiation of the delay clock.
  3. State-independent delay with state-dependent interruption: from the example of [3], consider a predator ($X$) prey ($Y$) model (including the events birth of prey, predation, and death of predators). Now assume that when prey give birth the juvenile prey requires time $\tau$ before becoming an adult. Also assume that juveniles are subject to predation, so they have a per-capita hazard of $\beta X$. Because the rate of predation upon juveniles is a function of state $X$ which may change before completion of the maturation phase $\tau$ due to other events, probability of interruption is not known at initiation.
  4. State-dependent uninterruptible delay: Consider a model where larvae compete for resources, so that the duration of time $\tau(X)$ required for maturation is explicitly a function of the state $X$, but that the larval stage is invulnerable to death, and so cannot be interrupted. An example is the damselfly model of [4] without death.
  5. State-dependent delay with state-independent interruption: Consider the damselfly model of [4]; a larvae's maturation delay $\tau(X)$ depends on the density of other larvae as they compete for a resource, but maturing larvae also suffer constant per-capita mortality $\mu$.
  6. State-dependent delay with state-dependent interruption: Consider complicating the damselfly model by making larval mortality state dependent, perhaps by introducing a predator class or cannibalism among the larvae.

For sampling, this means we want to associate a clock with two distributions.

## Sampling a Delayed Clock

Create a sampler to support delays by passing the `support_delayed=true` flag. For example:
```julia
builder = SamplerBuilder(Symbol, Float64;
                         support_delayed=true,
                         method=FirstToFireMethod(),
                         common_random=true)
sampler = SamplingContext(builder, rng)
```

Enable a delayed clock using a pair of distributions.

```julia
enable!(ctx, :recover, Exponential(0.5) => Dirac(5.0))
```

The first is the time until the reaction starts and the second is the time from the start to the end. Either or both can be of fixed duration or have distributions in time.

In the regular main loop, instead of calling `next(sampler)`, call `next_delayed(sampler)`.

```julia
when, which, phase = next_delayed(sampler)
```
The `phase` is one of three symbols.

 - `:regular` for clocks with no delays.
 - `:initiate` for the first phase of a delayed clock.
 - `:complete` for the second phase of a delayed clock.

You can `delete!(sampler, clock)` at any time, and it will delete either the initial and final, or just the final phase. It deletes whatever is currently enabled for that clock.

When the main loop fires, it should specify which clock it's firing by adding the phase.

```julia
fire!(ctx, which, phase, when) 
```

The initial clock isn't hidden from the user because some simulations need to be notified when the initial clock fires.

## Compatibility

 - Common random numbers and likelihoods are compatible with delayed clocks.
 - There isn't support for using arrays of distributions for vectorized importance sampling.
 - If you check for an enabled clock, you will see that its key includes the clock phase.
 