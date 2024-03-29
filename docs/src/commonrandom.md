# Common Random Numbers

## Introduction

If you set up the same model and run it with different initial random number generator states, then it will create a set of trajectories. Fleck sees these as a sequence of clock events and times of those events. You are usually interested in some summary outcomes of a simulation, such as the total time to a goal or the number of events. This summary outcome is a predictable function of the trajectories. We often want to ask how the goal function depends on simulation parameters, and that can be difficult to determine because each trajectory gives an individual value, and the set of trajectories gives an estimate that can have wide variance.

What we want is a variance reduction technique. Common random numbers are a variance reduction technique that enables you to use fewer simulation runs to compare the effect of different simulation parameters on the outcome. There are several other variance reduction techniques, such as antithetic variates and importance sampling, but let's look at common random numbers in Fleck.

## Using Common Random Numbers in Fleck

Fleck implements common random numbers by recording the state of the random number generator every time a clock is enabled. There are other ways to do this, but this one works with the [CombinedNextReaction](@ref) and [FirstToFire](@ref) samplers. The workflow you would use looks notionally like:

  1. Create a sampler.
  2. Wrap it in a [CommonRandomRecorder](@ref).
  3. Run a lot of simulations. Run `reset!(recorder)` after each simulation.
  4. For every parameter set to try, run it the same way, using `reset!` after each run.
  5. Compare outcomes.
