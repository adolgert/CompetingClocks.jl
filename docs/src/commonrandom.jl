# # Common Random Numbers

# ## Introduction

# If you set up the same model and run it with different initial random number generator (RNG) states, then it will create a set of trajectories. Fleck sees these as a sequence of clock events and times of those events. You are usually interested in some summary outcomes of a simulation, such as the total time to a goal or the number of events. This summary outcome is a predictable function of the trajectories. We often want to ask how the goal function depends on simulation parameters, and that can be difficult to determine because each trajectory gives an individual value, and the set of trajectories gives an estimate that can have wide variance.

# What we want is a [variance reduction](https://en.wikipedia.org/wiki/Variance_reduction) technique. Common random numbers (CRN) are a variance reduction technique that enables you to use fewer simulation runs to compare the effect of different simulation parameters on the outcome. There are several other variance reduction techniques, such as antithetic variates and importance sampling, but let's look at common random numbers in Fleck.

# If you estimate a value with $n$ independent trajectories, then the bias of the estimate is proportional to $1/\sqrt{n}$ in most cases. If you want to distinguish the effect of changing a parameter, then the estimate must be precise enough that you can see the difference. It is common to use millions of trajectories. On the other hand, CRN means that you can produce $n=100$ trajectories, with significant bias in the estimate, and still see the effect of changing a parameter.

# CRN works well when the sample path is similar from run to run. If two runs use completely different events, then there will be too little overlap. If the causal chain of which events affect other events changes, that can be a problem, too. In most cases, people try CRN and see if it helps.


# ## Using Common Random Numbers in Fleck

# Fleck implements common random numbers by recording the state of the random number generator every time a clock is enabled. There are other ways to do this, but this one works with the [CombinedNextReaction](@ref) and [FirstToFire](@ref) samplers. The workflow you would use looks notionally like:

#   1. Create a sampler.
#   2. Wrap it in a [CommonRandomRecorder](@ref).
#   3. Run a lot of simulations in order to explore and record all possible clock states. Run `reset!(recorder)` after each simulation.
#   4. For every parameter set to try, run it the same way, using `reset!` after each run.
#   5. Compare outcomes.

# Because the `CommonRandomRecorder` stores the state of the random number generator at each step, it works best with random number generators that have small state, such as Xoshiro on a linear congruential generator (LCG).

struct MakeModel end #hide
modify_model!(model, param_idx) = model #hide
run_simulation(model, sampler) = nothing #hide
using Random: Xoshiro
using Fleck
example_clock = (3, 7)  # We will use clock IDs that are a tuple of 2 integers.
model = MakeModel()
sampler = FirstToFire{typeof(example_clock),Float64}()
crn_sampler = CommonRandomRecorder(sampler, typeof(example_clock), Xoshiro)
for trial_idx in 1:100
    run_simulation(model, crn_sampler)
    reset!(crn_sampler)
end
for param_idx in 1:10
    each_model = modify_model!(model, param_idx)
    run_simulation(each_model, crn_sampler)
    reset!(crn_sampler)
end

# ## Multithreading

# A joy of using simulations is how easy it is to parallelize simulation runs across tasks. That can be a challenge for the `CommonRandomRecorder` because it continues to observe and record new RNG states as it comes across them. That will result in divergence between behavior on different threads. For that reason, it is possible to [freeze](@ref) a `CommonRandomRecorder`. It will stop recording states, so make sure to first prime it with lots of simulation runs, and then freeze the recorder and use that as the sampler for multiple simulations on multiple threads.

using Random: Xoshiro
using Fleck
example_clock = (3, 7)  # We will use clock IDs that are a tuple of 2 integers.
model = MakeModel()
sampler = FirstToFire{typeof(example_clock),Float64}()
crn_sampler = CommonRandomRecorder(sampler, typeof(example_clock), Xoshiro)
for trial_idx in 1:100
    run_simulation(model, crn_sampler)
    reset!(crn_sampler)
end
for thread_idx in 1:10
    frozen_crn = freeze(crn_sampler)
    ## start a simulation run on this thread with frozen_crn.
end


# ## Checking effectiveness of Common Random Numbers

# If your simulation has a large sample space, CRN may not help. We run a first set of simulations in order to record the state of the system for lots of different clocks and different multiplicities of clock events. If that worked well, then subsequent runs of the simulation will re-use draws from the random number generator. If there are a lot of events which are needed but haven't been recorded, those misses are a sign that CRN is unlikely to reduce variance much for this simulation.

# We check this by checking the [misscount](@ref) during later runs of the simulation under CRN. If that miss count is high, we can look into which clocks are firing that didn't previously fire by iterating over the [misses](@ref), which are pairs of (clock key, number of misses for that clock).

# The final word on effectiveness of CRN is to look at the variance of summary outcomes for runs with and without CRN. The CRN will, in general, slow down a sampler, but it should mean that many fewer runs are required to distinguish the effect of changes in system parameters.
