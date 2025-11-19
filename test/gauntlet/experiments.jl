using .TravelModel

"""
    mark_calibration_nonparametric_bootstrap(commands, sampler_cnt, rng)

This resamples with replacement in order to mimic repeated sampling from the
same data-generating process.
"""
function mark_calibration_nonparametric_bootstrap(commands, sampler_cnt, rng)
    B = 1000  # more stable than 100
    boot = zeros(Float64, B)

    # Build one “dataset” of draws once
    samplers, final_time = parallel_replay(commands, sampler_cnt, rng)
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


function mark_calibration_single()
    rng = [Xoshiro(98327423 + 298432*i) for i in 1:Threads.maxthreadid()]
    model = Travel(5, TravelGraph.complete, TravelMemory.forget, rng[1])
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng[1])

    sampler_cnt = 10000
    mark_calibration_nonparametric_bootstrap(commands, sampler_cnt, rng)

    samplers, final_time = parallel_replay(commands, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)
    mark_score = mark_calibration_conditional_time(draws, distributions)

    pvalue = mark_calibration_permutation_test(draws, distributions, mark_score, rng[1])
    @show mark_score
    @show pvalue
end
