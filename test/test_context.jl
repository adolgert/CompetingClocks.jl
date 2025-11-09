# For this macro, need to import it into main for it to work.
import InteractiveUtils: @code_typed



@testset "SamplerContext doesn't penalize missing pieces" begin
    using CompetingClocks
    using Random
    using Distributions

    K = Int64
    T = Float64
    rng = Xoshiro(899987987)
    builder = SamplerBuilder(K, T; method=FirstToFireMethod())
    sampler = SamplingContext(builder, rng)
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity))
    end
    # Assert that when compiled, the if-then statements in the context
    # are compiled away!
    res = @code_typed enable!(sampler, 1, Exponential(0.5), 0.0)
    ci = first(res)
    branch_count = count(expr -> isa(expr, Core.GotoIfNot), ci.code)
    @test branch_count == 0
end


@safetestset context_life_cycle = "Context works with everybody's life cycle" begin
    using CompetingClocks
    using Random
    using Distributions

    K = Int64
    T = Float64
    for SamplerType in [FirstToFireMethod(), DirectMethod(), FirstReactionMethod(), NextReactionMethod()]
        rng = Xoshiro(90422342)
        builder = SamplerBuilder(K, T; method=SamplerType)
        context = SamplingContext(builder, rng)
        test_enable = Set{Int64}()
        for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
            enable!(context, clock_id, Exponential(propensity))
            push!(test_enable, clock_id)
        end
        when, which = next(context)
        fire!(context, which, when)
        delete!(test_enable, which)
        when, which = next(context)
        @test when > 0.0
        @test 1 <= which
        @test which <= 5
        @test which ∈ test_enable
        @test enabled(context) == test_enable
        for try_one in test_enable
            @test isenabled(context, try_one)
        end
    end
end
