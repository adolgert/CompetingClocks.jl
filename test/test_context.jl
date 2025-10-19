# For this macro, need to import it into main for it to work.
import InteractiveUtils: @code_typed


@safetestset samplerbuilder_keyclassifier = "SamplerBuilder key classifier" begin
    using CompetingClocks: make_key_classifier
    using Distributions
    testset = Dict(
        :high => (k, d) -> k[1] == :mountain,
        :low => (k, d) -> k[1] == :valley
    )
    func = make_key_classifier(testset)
    @test func((:mountain, 17), Exponential()) == :high
    @test func((:valley, "saguaro"), Weibull()) == :low
    @test_throws ErrorException func((:forest, 32), Exponential())
end


@safetestset samplerbuilder_construct_single = "SamplerBuilder single construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple,Float64)
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    @test isempty(scbuild.group)
    add_group!(scbuild, :capistrano => (k,d) -> true; sampler_spec=:firsttofire)
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_once = "SamplerBuilder single construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple,Float64; sampler_spec=:firsttofire)
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    @test !isempty(scbuild.group)
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_multi = "SamplerBuilder multi construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple,Float64)
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    add_group!(scbuild, :sparky => (x,d) -> x[1] == :recover, sampler_spec=(:nextreaction,))
    add_group!(scbuild, :forthright=>(x,d) -> x[1] == :infect)
    scresult = build_sampler(scbuild)
    @test scresult isa MultiSampler
end


@testset "SamplerContext doesn't penalize missing pieces" begin
    using CompetingClocks
    using Random
    using Distributions

    K = Int64
    T = Float64
    rng = Xoshiro(899987987)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire)
    sampler = SamplingContext(builder, rng)
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0)
    end
    # Assert that when compiled, the if-then statements in the context
    # are compiled away!
    res = @code_typed enable!(sampler, 1, Exponential(0.5), 0.0, 0.0)
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
    for SamplerType in [:firsttofire, :direct, :firstreaction, :nextreaction]
        rng = Xoshiro(90422342)
        builder = SamplerBuilder(K, T; sampler_spec=SamplerType)
        context = SamplingContext(builder, rng)
        test_enable = Set{Int64}()
        for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
            enable!(context, clock_id, Exponential(propensity), 0.0, 0.0)
            push!(test_enable, clock_id)
        end
        when, which = next(context, 0.0)
        fire!(context, which, when)
        delete!(test_enable, which)
        when, which = next(context, when)
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
