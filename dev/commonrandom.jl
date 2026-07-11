# # Common Random Numbers

# ## Introduction

# If you set up the same model and run it with different initial random number generator (RNG) states, then it will create a set of trajectories. CompetingClocks sees these as a sequence of clock events and times of those events. You are usually interested in some summary outcomes of a simulation, such as the total time to a goal or the number of events. This summary outcome is a predictable function of the trajectories. We often want to ask how the goal function depends on simulation parameters, and that can be difficult to determine because each trajectory gives an individual value, and the set of trajectories gives an estimate that can have wide variance.

# What we want is a [variance reduction](https://en.wikipedia.org/wiki/Variance_reduction) technique. Common random numbers (CRN) are a variance reduction technique that enables you to use fewer simulation runs to compare the effect of different simulation parameters on the outcome. There are several other variance reduction techniques, such as antithetic variates and importance sampling, but let's look at common random numbers in CompetingClocks.

# CRN works well when the sample path is similar from run to run. If two runs use completely different events, then there will be too little overlap. If the causal chain of which events affect other events changes, that can be a problem, too. In most cases, people try CRN and see if it helps.


# ## Using Common Random Numbers in CompetingClocks

# Each sampler OWNS its randomness as a family of per-clock keyed streams (see
# [CompetingClocks.KeyedStreams](@ref)). A draw belongs to a CLOCK — and to that
# clock's occurrence count — not to a position in a global call sequence. The
# consequence is the whole of common random numbers: two samplers built from the
# **same seed** draw the identical firing time for the same clock, even when their
# events fire in a different order, because each clock's stream is seeded once
# from the clock and advances only when that clock draws.

# There is therefore no recorder to prime, no replay mode, and no miss count. To
# compare a baseline against a perturbed parameter with common random numbers, you
# build both from the same seed and run them:

struct MakeModel end #hide
modify_model!(model, param_idx) = model #hide
run_simulation(model, sampler) = 0.0 #hide
using Random: Xoshiro
using CompetingClocks
example_clock = (3, 7)  # We will use clock IDs that are a tuple of 2 integers.
model = MakeModel()
(Key, Time) = (typeof(example_clock), Float64)

# The context's rng is used only to choose the sampler's stream seed, so passing
# the SAME rng seed to two contexts couples them. Compare a baseline outcome to a
# perturbed one, reusing the same seed across the paired runs:
baseline = Float64[]
perturbed = Float64[]
for trial_idx in 1:100
    ## The same seed drives both runs of the pair, which couples them.
    base_ctx = SamplingContext(SamplerBuilder(Key, Time), Xoshiro(trial_idx))
    pert_ctx = SamplingContext(SamplerBuilder(Key, Time), Xoshiro(trial_idx))
    push!(baseline, run_simulation(model, base_ctx))
    push!(perturbed, run_simulation(modify_model!(model, 1), pert_ctx))
end

# The paired difference `perturbed .- baseline` has far smaller variance than if
# the two runs had used independent seeds, because each clock's draw is shared
# between the pair. That is the variance reduction, and it falls out of stream
# ownership rather than any bookkeeping.

# ## Cloning and re-keying

# Two lower-level primitives express the same idea inside one run:
#
#   * `clone(sampler)` makes a full-state copy — clock state AND the keyed streams
#     (generator states and occurrence counts). The clone and the original race
#     the identical clocks to the identical times: they are coupled.
#   * `rekey_streams!(sampler, seed)` re-seeds the streams, forgetting every live
#     generator, so a clone can be DEcoupled to explore a divergent future.
#
# Splitting a run into independent continuations is exactly a `clone` followed by
# a `rekey_streams!` on each copy. The [Randomness Ownership](randomness.md) page
# explains the stream mechanism, `split!`, and the low-level coupling recipe in
# detail.

# ## Checking effectiveness of Common Random Numbers

# If your simulation has a large sample space, CRN may not help: when the two runs
# take genuinely different events, few clocks share a coupled draw. The final word
# on effectiveness is to look at the variance of summary outcomes for the paired
# (same-seed) runs versus independent-seed runs. If the paired variance is much
# smaller, CRN is helping and you need far fewer runs to distinguish the effect of
# a parameter change.
