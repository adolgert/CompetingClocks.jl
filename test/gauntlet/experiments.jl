using .TravelModel
using UnitTestDesign


struct SamplerSUT
    travel::TravelConfig
    state_cnt::Int
    step_cnt::Int
end


function rng_set(single_rng)
    rng = Vector{Xoshiro}(undef, Threads.maxthreadid())
    rng[1] = single_rng
    for rng_gen in 2:length(rng)
        rng[rng_gen] = jump!(rng[rng_gen - 1])
    end
    return rng
end


function mark_calibration_single(smethod::SamplerSpec, sut::SamplerSUT, rng::Vector)
    sampler = smethod(Int,Float64)
    commands = travel_commands(sut.step_cnt, sut.state_cnt, sut.travel, rng[1])

    sampler_cnt = 10000
    mark_calibration_nonparametric_bootstrap(commands, sampler_cnt, rng)

    samplers, final_time = parallel_replay(commands, sampler, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)
    mark_score = mark_calibration_conditional_time(draws, distributions)

    pvalue = mark_calibration_permutation_test(draws, distributions, mark_score, rng[1])
    return (; test="mark-calibration", mark_score, pvalue)
end


function doob_meyer_single(smethod::SamplerSpec, sut::SamplerSUT, rng::Vector)
    sampler = smethod(Int,Float64)
    commands = travel_commands(sut.step_cnt, sut.state_cnt, sut.travel, rng[1])
    distributions = final_enabled_distributions(commands)

    sampler_cnt = 10000
    draws, final_time = retrieve_draws(commands, sampler, sampler_cnt, rng)
    dm_test = doob_meyer(draws, collect(values(distributions)), final_time)
    return dm_test
end


function two_sample_ad_single(smethod::SamplerSpec, sut::SamplerSUT, rng::Vector)
    sampler = smethod(Int,Float64)
    commands = travel_commands(sut.step_cnt, sut.state_cnt, sut.travel, rng[1])
    distributions = final_enabled_distributions(commands)

    sampler_cnt = 10000
    compare = FirstReactionMethod()(Int, Float64)
    draws_a, final_time = retrieve_draws(commands, compare, sampler_cnt, rng)
    draws_b, final_time = retrieve_draws(commands, sampler, sampler_cnt, rng)

    return ad_two_sample(draws_a, draws_b, collect(keys(distributions)); verbose=true)
end


function experiment_set(sampler_method, sut, rng_single)
    rng = rng_set(rng_single)
    marks = mark_calibration_single(sampler_method, sut, rng)
    times = doob_meyer_single(sampler_method, sut, rng)
    comps = two_sample_ad_single(sampler_method, sut, rng)
    push!(comps, marks)
    push!(comps, times)
    return comps
end


function experiment_range()
    memory = [TravelMemory.forget,]
    graph = [TravelGraph.cycle, TravelGraph.complete,]
    dist = [TravelRateDist.exponential, TravelRateDist.general]
    count = [TravelRateCount.destination, TravelRateCount.pair]
    delay = [TravelRateDelay.none, TravelRateDelay.right]
    sampler_spec = [FirstReactionMethod(), FirstToFireMethod()]
    step_cnt = [5, 6]
    state_cnt = [4, 5]
    configurations = full_factorial(
        memory, graph, dist, count, delay, collect(1:length(sampler_spec)), step_cnt, state_cnt
        )
    arrangements = Vector{Tuple{SSA,SamplerSUT}}(undef, length(configurations))
    for idx in eachindex(configurations)
        configuration = configurations[idx]
        config = TravelConfig(configuration[1:5]...)
        sampler = sampler_spec[configuration[6]]
        step_cnt = configuration[7]
        state_cnt = configuration[8]
        arrangments[idx] = (sampler, SamplerSUT(config, state_cnt, step_cnt))
    end
    return arrangements
end


function run_experiments()
    rng_single = Xoshiro(882342987)
    configurations = experiment_range()
    results = Vector{Any}(undef, length(configurations))
    for gen_idx in eachindex(configurations)
        sampler, sut = configurations[gen_idx]
        results[gen_idx] = experiment_set(sampler, sut, rng_single)
    end
    scores = Vector{Tuple{Float64,Int}}(undef, length(results))
    for score_idx in eachindex(results)
        scores[score_idx] = (minimum(x.pvalue for x in results[score_idx]), score_idx)
    end
    sort!(scores)
    println("lowest scores")
    println(scores[begin:begin + 10])
end
