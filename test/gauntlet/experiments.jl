using .TravelModel


function rng_set(single_rng)
    rng = Vector{Xoshiro}(undef, Threads.maxthreadid())
    rng[1] = single_rng
    for rng_gen in 2:length(rng)
        rng[rng_gen] = jump!(rng[rng_gen - 1])
    end
    return rng
end


function mark_calibration_single(single_rng)
    rng = rng_set(single_rng)
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    step_cnt = 5
    state_cnt = 5
    commands = travel_commands(step_cnt, state_cnt, config, rng[1])

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


function doob_meyer_single(single_rng)
    rng = rng_set(single_rng)
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    step_cnt = 5
    state_cnt = 5
    commands = travel_commands(step_cnt, state_cnt, config, rng[1])
    distributions = final_enabled_distributions(commands)

    sampler_cnt = 10000
    draws, final_time = retrieve_draws(commands, sampler_cnt, rng)
    dm_test = doob_meyer(draws, collect(values(distributions)), final_time)
    @show dm_test.test
    @show dm_test.pvalue
end


function two_sample_ad_single(single_rng)
    rng = rng_set(single_rng)
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    step_cnt = 5
    state_cnt = 5
    commands = travel_commands(step_cnt, state_cnt, config, rng[1])
    distributions = final_enabled_distributions(commands)

    sampler_cnt = 10000
    draws_a, final_time = retrieve_draws(commands, sampler_cnt, rng)
    draws_b, final_time = retrieve_draws(commands, sampler_cnt, rng)

    ad_two_sample(draws_a, draws_b, collect(keys(distributions)); verbose=true)
end
