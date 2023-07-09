# This file tests the performance of the CombinedNextReaction type.
# It compares the performance of that sampler against the performance
# of the classic NextReaction and ModifiedNextReaction types. We are
# guarding against errors that come from mishandling types in Julia.
using Distributions
using Fleck
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
    pure_time = 0
    dual_time = 0
    # We run a few iterations in order to account for compilation time.
    for burn_one in 1:3
        pure_nr = NextReaction{Int}()
        dual_nr = CombinedNextReaction{Int}()
        pure_rng = Xoshiro(342432)
        pure_time = @timed sample_a_while(pure_nr, common_distribution, pure_rng)
        dual_rng = Xoshiro(342432)
        dual_time = @timed sample_a_while(dual_nr, common_distribution, dual_rng)
    end
    return (pure_time.time, dual_time.time)
end


function compare_with_modified_next_reaction()
    common_distribution = Exponential()
    pure_time = 0
    dual_time = 0
    # We run a few iterations in order to account for compilation time.
    for burn_one in 1:3
        pure_nr = ModifiedNextReaction{Int}()
        dual_nr = CombinedNextReaction{Int}()
        pure_rng = Xoshiro(342432)
        pure_time = @timed sample_a_while(pure_nr, common_distribution, pure_rng)
        dual_rng = Xoshiro(342432)
        dual_time = @timed sample_a_while(dual_nr, common_distribution, dual_rng)
    end
    return (pure_time.time, dual_time.time)
end


println(compare_with_next_reaction())
println(compare_with_modified_next_reaction())
