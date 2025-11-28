using Base.Threads: @threads, maxthreadid, threadid
using CompetingClocks: hazard

"""
    mark_calibration_brier(distributions, fired_clock, when)

Calculating one component of the Brier score. The whole Brier score is

```math
    \\frac{1}{N}\\sum_i^N \\sum_i^R (f_{ti}-o_{ti})^2
```

where ``o_{ti}`` is 1 for the fired clock and 0 otherwise. Here ``f_{ti}`` is
the hazard of one clock divided by the total hazard. This function calculates
the above for one outcome ``i``.
"""
function mark_calibration_brier(distributions, fired_clock, when)
    total_hazard = 0.0
    hazards = Dict{Int,Float64}()
    for (clock, ds) in distributions
        h = hazard(ds.d, ds.enabling_time, when)
        hazards[clock] = h
        total_hazard += h
    end
    total_hazard == 0.0 && return 0.0

    s = 0.0
    for (clock, h) in hazards
        p = h / total_hazard
        y = clock == fired_clock ? 1.0 : 0.0
        s += (p - y)^2
    end
    # Divide by distributions because this it the number of classes
    # for multiclass Brier.
    return s / length(distributions)
end


"""
    mark_calibration_conditional_time(draws, distributions)

This tests whether ``P[E|T]`` is correct. It's the probability that the right
clock fires given when a clock fired. It's a simple value because it's the hazard
for the clock that fired divided by the sum of all hazards at the time that clock
fired.
"""
function mark_calibration_conditional_time(draws, distributions)
    # Use thread-local storage for reduction
    partial_sums = zeros(Float64, Threads.maxthreadid())

    @threads for i in eachindex(draws)
        tid = Threads.threadid()
        partial_sums[tid] += mark_calibration_brier(distributions, draws[i][1], draws[i][2])
    end

    return sum(partial_sums) / length(draws)
end


"""
    mark_calibration_nonparametric_bootstrap(commands, sampler_cnt, rng)

This resamples with replacement in order to mimic repeated sampling from the
same data-generating process.
"""
function mark_calibration_nonparametric_bootstrap(commands, sampler, sampler_cnt, rng)
    B = 1000  # more stable than 100
    boot = zeros(Float64, B)

    # Build one “dataset” of draws once
    samplers, final_time = parallel_replay(commands, sampler, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)

    for b in 1:B
        idx = rand(rng[1], 1:length(draws), length(draws))   # sample indices w/ replacement
        bs_draws = [draws[i] for i in idx]
        boot[b] = mark_calibration_conditional_time(bs_draws, distributions)
    end

    sort!(boot)
    lo = boot[ceil(Int, 0.025*(B+1))]
    hi = boot[floor(Int, 0.975*(B+1))]
    return @show (lo=lo, hi=hi)
end


"""
    mark_calibration_resimulation(commands, sampler_cnt, rng)

This generates a lot of Brier scores in order to see variation in what
the sampler produces.
"""
function mark_calibration_resimulation(commands, sampler_cnt, rng)
    bootstrap_cnt = 1000
    bootstraps = zeros(Float64, bootstrap_cnt)
    for boot_idx in 1:bootstrap_cnt
        samplers, final_time = parallel_replay(commands, sampler_cnt, rng)
        draws = sample_samplers(samplers, final_time, rng)
        distributions = final_enabled_distributions(commands)
        mark_score = mark_calibration_conditional_time(draws, distributions)
        bootstraps[boot_idx] = mark_score
    end
    sort!(bootstraps)
    lo = bootstraps[begin]
    hi = bootstraps[end]
    edge = bootstraps[Int(round(bootstrap_cnt * 0.025))]
    return @show (lo, edge, hi)
end


"""
A permutation test of the Brier score.
Computes a one-sided p-value.

```math
p = \\frac{1+\\mbox{count}(b \\le obs)}{B+1}
```
Lower Brier is better.
"""
function mark_calibration_permutation_test(draws, distributions, obs, rng::AbstractRNG)
    B = 2000
    perm = zeros(Float64, B)
    for b in 1:B
        d = copy(draws)
        jumble!(d, rng)  # break the pairing
        perm[b] = mark_calibration_conditional_time(d, distributions)
    end
    pval = (1 + count(x -> x <= obs, perm)) / (B + 1)
    return pval
end
