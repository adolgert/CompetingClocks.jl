# This file tests the performance of the CombinedNextReaction type.
# It compares the performance of that sampler against the performance
# of the classic NextReaction and ModifiedNextReaction types. We are
# guarding against errors that come from mishandling types in Julia.
using Distributions
using CompetingClocks
using Random


function sample_a_while(sampler, distribution, rng)
    for create_idx in 1:100
        enable!(sampler, create_idx, distribution, 0.0, 0.0, rng)
    end
    when = 0.0
    for select_idx in 1:10000
        when, chosen = next(sampler, when, rng)
        disable!(sampler, chosen, when)
        enable!(sampler, chosen, distribution, when, when, rng)
    end
    return when
end


function compare_with_next_reaction()
    common_distribution = InverseGamma()
    pure_result = nothing
    dual_result = nothing
    # We run a few iterations in order to account for compilation time.
    iter_cnt = 0
    while pure_result === nothing || dual_result === nothing
        pure_nr = NextReaction{Int}()
        dual_nr = CombinedNextReaction{Int,Float64}()
        pure_rng = Xoshiro(342432)
        dual_rng = Xoshiro(342432)
        pure_time = @timed sample_a_while(pure_nr, common_distribution, pure_rng)
        dual_time = @timed sample_a_while(dual_nr, common_distribution, dual_rng)
        iter_cnt += 1
        if iter_cnt > 1
            if pure_time.gctime == 0.0
                pure_result = pure_time
            end
            if dual_time.gctime == 0.0
                dual_result = dual_time
            end
        end
    end
    return (pure_result, dual_result)
end


function compare_with_modified_next_reaction()
    common_distribution = Exponential()
    pure_time = @timed sin(0.1)
    dual_time = @timed cos(0.1)
    # We run a few iterations in order to account for compilation time.
    for burn_one in 1:4
        pure_nr = ModifiedNextReaction{Int}()
        dual_nr = CombinedNextReaction{Int,Float64}()
        pure_rng = Xoshiro(342432)
        dual_rng = Xoshiro(342432)
        pure_time = @timed sample_a_while(pure_nr, common_distribution, pure_rng)
        dual_time = @timed sample_a_while(dual_nr, common_distribution, dual_rng)
    end
    return (pure_time, dual_time)
end


pure, dual = compare_with_next_reaction()
println("next reaction")
println(pure)
println("combined next reaction")
println(dual)
mnrm, combo = compare_with_modified_next_reaction()
println("modified next reaction")
println(mnrm)
println("combined next reaction")
println(combo)
