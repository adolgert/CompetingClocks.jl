using .TravelModel

function single_mark_calibration()
    rng = Xoshiro(98327423)
    model = Travel(2, TravelGraph.path, TravelMemory.forget, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng)

    sampler_cnt = 10000
    samplers, final_time = parallel_replay(commands, sampler_cnt, 2394723)
    draws = sample_samplers(samplers, final_time, rng)
    distributions = final_enabled_distributions(commands)
    mark_score = mark_calibration_conditional_time(draws, distributions)
    @show mark_score
end
