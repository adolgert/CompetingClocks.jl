using .TravelModel
using UnitTestDesign


function rng_set(single_rng)
    rng = Vector{Xoshiro}(undef, Threads.maxthreadid())
    rng[1] = single_rng
    for rng_gen in 2:length(rng)
        rng[rng_gen] = jump!(rng[rng_gen - 1])
    end
    return rng
end


function mark_calibration_single(smethod, sut::SamplerSUT, rng::Vector)
    sampler = smethod(Int,Float64)
    commands = travel_commands(sut.step_cnt, sut.state_cnt, sut.travel, rng[1])

    sampler_cnt = 10000
    mark_calibration_nonparametric_bootstrap(commands, sampler, sampler_cnt, rng)

    samplers, final_time = parallel_replay(commands, sampler, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)
    mark_score = mark_calibration_conditional_time(draws, distributions)

    pvalue = mark_calibration_permutation_test(draws, distributions, mark_score, rng[1])
    return (; test="mark-calibration", mark_score, pvalue)
end


function doob_meyer_single(smethod, sut::SamplerSUT, rng::Vector)
    sampler = smethod(Int,Float64)
    commands = travel_commands(sut.step_cnt, sut.state_cnt, sut.travel, rng[1])
    distributions = final_enabled_distributions(commands)

    sampler_cnt = 10000
    draws, final_time = retrieve_draws(commands, sampler, sampler_cnt, rng)
    dm_test = doob_meyer(draws, collect(values(distributions)), final_time)
    return dm_test
end


function two_sample_ad_single(smethod, sut::SamplerSUT, rng::Vector)
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


function experiment_focused()
    arrangements = Vector{Tuple{CompetingClocks.SamplerSpec,SamplerSUT}}()
    push!(arrangements,
        (
            RejectionMethod(),
            SamplerSUT(
                TravelConfig(
                    TravelMemory.forget,
                    TravelGraph.cycle,
                    TravelRateDist.exponential,
                    TravelRateCount.destination,
                    TravelRateDelay.none
                ),
                5,
                10_000
            )
        )
    )
    return arrangements
end


function experiment_exponential(faster::Bool)
    graph = [TravelGraph.cycle, TravelGraph.complete,]
    count = [TravelRateCount.destination, TravelRateCount.pair,]
    sampler_spec = [
        RejectionMethod(), PartialPropensityMethod(), DirectMethod(:keep, :tree),
        DirectMethod(:keep, :array), DirectMethod(:remove, :tree),
        DirectMethod(:remove, :array),
        ]
    state_cnt = [4, 5, 6]
    design = faster ? all_pairs : full_factorial
    configurations = design(
        graph, count, collect(1:length(sampler_spec)),
        state_cnt
        )
    arrangements = Vector{Tuple{CompetingClocks.SamplerSpec,SamplerSUT}}(undef, length(configurations))
    memory = TravelMemory.forget
    dist = TravelRateDist.exponential
    delay = TravelRateDelay.none
    for idx in eachindex(configurations)
        configuration = configurations[idx]
        config = TravelConfig(memory, configuration[1], dist, configuration[2], delay)
        sampler = sampler_spec[configuration[3]]
        state_cnt = configuration[4]
        arrangements[idx] = (sampler, SamplerSUT(config, state_cnt, 10000))
    end
    return arrangements
end

function experiment_range(faster::Bool)
    memory = [TravelMemory.forget, TravelMemory.remember,]
    graph = [TravelGraph.cycle, TravelGraph.complete,]
    dist = [TravelRateDist.exponential, TravelRateDist.general,]
    count = [TravelRateCount.destination, TravelRateCount.pair,]
    delay = [TravelRateDelay.none, TravelRateDelay.right,]
    sampler_spec = [
        FirstReactionMethod(), FirstToFireMethod(), NextReactionMethod(),
        ]
    state_cnt = [4, 5, 6]
    design = faster ? all_pairs : full_factorial
    configurations = design(
        memory, graph, dist, count, delay, collect(1:length(sampler_spec)),
        state_cnt
        )
    arrangements = Vector{Tuple{CompetingClocks.SamplerSpec,SamplerSUT}}(undef, length(configurations))
    for idx in eachindex(configurations)
        configuration = configurations[idx]
        config = TravelConfig(configuration[1:5]...)
        sampler = sampler_spec[configuration[6]]
        state_cnt = configuration[7]
        arrangements[idx] = (sampler, SamplerSUT(config, state_cnt, 10_000))
    end
    return arrangements
end


function run_experiments(faster=false)
    echo = !faster
    rng_single = Xoshiro(882342987)
    configurations1 = experiment_range(faster)
    configurations2 = experiment_exponential(faster)
    configurations = vcat(configurations2, configurations1)
    # configurations = experiment_focused()
    echo && println("There are $(length(configurations)) configurations.")
    results = Vector{Any}(undef, length(configurations))
    for gen_idx in eachindex(configurations)
        sampler_spec, sut = configurations[gen_idx]
        echo && println("spec $sampler_spec sut $sut")
        results[gen_idx] = collect_data_single(sampler_spec, sut, rng_single)
    end
    scores = Vector{Tuple{Float64,Int}}(undef, length(results))
    for score_idx in eachindex(results)
        mmm = minimum(x.pvalue for x in results[score_idx])
        # We take the minimum value across a set of clocks, so it is
        # a minimum of a different number of values. Each of the values
        # should be uniformly distributed, so the min is a convolution.
        # Using this adjusts each min for the number of metrics.
        adjusted = 1 - (1 - mmm)^length(results[score_idx])
        scores[score_idx] = (adjusted, score_idx)
    end
    sort!(scores)
    echo && println("=" ^ 80)
    echo && println("lowest scores")
    echo && println("=" ^ 80)
    succeed = true
    for examine in 1:5
        value, config_idx = scores[examine]
        config = configurations[config_idx]
        echo && println("value $value")
        echo && println("config $config")
        res_metrics = results[config_idx]
        group = Dict{Tuple{String,Int},Vector{Float64}}()
        for res in res_metrics
            echo && println("metric $(res.name) $(res.pvalue) $(res.clock) $(res.count)")
            group[(res.name, res.clock)] = [res.pvalue]
        end
        echo && println("=" ^ 80)
        sampler_spec, sut = configurations[config_idx]
        for i in 1:5
            rep_metrics = collect_data_single(sampler_spec, sut, rng_single)
            echo && println("-"^80)
            for res in rep_metrics
                echo && println("metric $(res.name) $(res.pvalue) $(res.clock) $(res.count)")
                push!(group[(res.name, res.clock)], res.pvalue)
            end
        end
        for (metid, metvals) in group
            m = median(metvals)
            @assert m > 0.05
            if m < 0.05
                succeed = false
            end
        end
        echo && println("=" ^ 80)
    end
    return succeed
end
