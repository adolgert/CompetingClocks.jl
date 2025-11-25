using .TravelModel


function mark_calibration_single()
    rng = [Xoshiro(98327423 + 298432*i) for i in 1:Threads.maxthreadid()]
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    model = Travel(5, config, rng[1])
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


function doob_meyer_single()
    rng = [Xoshiro(98327423 + 298432*i) for i in 1:Threads.maxthreadid()]
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    model = Travel(5, config, rng[1])
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng[1])

    sampler_cnt = 10000
    samplers, final_time = parallel_replay(commands, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)
    forget_clock = [x[2] for x in draws]
    dm_test = doob_meyer(forget_clock, collect(values(distributions)), final_time)
    @show dm_test
    @show pvalue(dm_test)
end


function two_sample_ad_single()
    rng = [Xoshiro(98327423 + 298432*i) for i in 1:Threads.maxthreadid()]
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.complete, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    model = Travel(5, config, rng[1])
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng[1])

    sampler_cnt = 10000
    draws_a = retrieve_draws(commands, sampler_cnt, rng)
    draws_b = retrieve_draws(commands, sampler_cnt, rng)

    distributions = final_enabled_distributions(commands)
    ad_two_sample(draws_a, draws_b, collect(keys(distributions)); verbose=true)
end
