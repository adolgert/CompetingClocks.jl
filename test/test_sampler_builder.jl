using SafeTestsets


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
    scbuild = SamplerBuilder(Tuple, Float64)
    @test isempty(scbuild.group)
    add_group!(scbuild, :capistrano => (k, d) -> true; method=FirstToFireMethod())
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_once = "SamplerBuilder single construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple, Float64; method=FirstToFireMethod())
    @test !isempty(scbuild.group)
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_multi = "SamplerBuilder multi construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple, Float64)
    add_group!(scbuild, :sparky => (x, d) -> x[1] == :recover; method=NextReactionMethod())
    add_group!(scbuild, :forthright => (x, d) -> x[1] == :infect)
    scresult = build_sampler(scbuild)
    @test scresult isa MultiSampler
end


@safetestset samplerbuilder_has_traits = "SamplerBuilder has_steploglikelihood/has_pathloglikelihood" begin
    using CompetingClocks: has_steploglikelihood, has_pathloglikelihood
    using CompetingClocks: CombinedNextReaction, DirectCall, MultipleDirect
    using CompetingClocks: EnabledWatcher, FirstToFire

    # has_steploglikelihood
    @test has_steploglikelihood(CombinedNextReaction) == true
    @test has_steploglikelihood(DirectCall) == true
    @test has_steploglikelihood(MultipleDirect) == true
    @test has_steploglikelihood(EnabledWatcher) == true
    @test has_steploglikelihood(FirstToFire) == false
    @test has_steploglikelihood(String) == false  # arbitrary type

    # has_pathloglikelihood
    @test has_pathloglikelihood(DirectCall) == true
    @test has_pathloglikelihood(MultipleDirect) == true
    @test has_pathloglikelihood(FirstToFire) == false
    @test has_pathloglikelihood(String) == false
end


@safetestset samplerbuilder_no_groups = "SamplerBuilder no groups" begin
    using CompetingClocks: SamplerBuilder, build_sampler, FirstToFire

    # When no groups are added, should default to FirstToFire
    builder = SamplerBuilder(Int64, Float64)
    @test isempty(builder.group)
    sampler = build_sampler(builder)
    @test sampler isa FirstToFire
end


@safetestset samplerbuilder_auto_select_path = "SamplerBuilder auto_select with path_likelihood" begin
    using CompetingClocks: SamplerBuilder, add_group!, build_sampler, DirectCall

    # With path_likelihood=true, should auto-select DirectMethod
    builder = SamplerBuilder(Int64, Float64; path_likelihood=true)
    add_group!(builder, :all => (k, d) -> true)  # no method specified
    sampler = build_sampler(builder)
    @test sampler isa DirectCall
end


@safetestset samplerbuilder_auto_select_step = "SamplerBuilder auto_select with step_likelihood" begin
    using CompetingClocks: SamplerBuilder, add_group!, build_sampler, CombinedNextReaction

    # With step_likelihood=true (but not path), should auto-select NextReactionMethod
    builder = SamplerBuilder(Int64, Float64; step_likelihood=true)
    add_group!(builder, :all => (k, d) -> true)  # no method specified
    sampler = build_sampler(builder)
    @test sampler isa CombinedNextReaction
end


@safetestset samplerbuilder_selector_error = "SamplerBuilder selector inconsistency error" begin
    using CompetingClocks: SamplerBuilder, add_group!

    # Adding a group without selector after one with selector should error
    builder = SamplerBuilder(Int64, Float64)
    add_group!(builder, :first => (k, d) -> true)
    @test_throws ErrorException add_group!(builder)  # no selector

    # Adding a group with selector after one without should also error
    builder2 = SamplerBuilder(Int64, Float64)
    add_group!(builder2)  # no selector (selector=nothing)
    @test_throws ErrorException add_group!(builder2, :second => (k, d) -> true)
end


@safetestset samplerbuilder_from_inclusion = "SamplerBuilder FromInclusion choose_sampler" begin
    using CompetingClocks: SamplerBuilder, add_group!, build_sampler
    using CompetingClocks: enable!, next
    using Distributions: Exponential
    using Random: Xoshiro

    rng = Xoshiro(123456)

    # Create a multi-sampler to exercise FromInclusion.choose_sampler
    builder = SamplerBuilder(Symbol, Float64)
    add_group!(builder, :fast => (k, d) -> k == :quick)
    add_group!(builder, :slow => (k, d) -> k == :gradual)
    sampler = build_sampler(builder)

    # Enable clocks - this exercises choose_sampler internally
    enable!(sampler, :quick, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, :gradual, Exponential(2.0), 0.0, 0.0, rng)

    @test length(sampler) == 2
    when, which = next(sampler, 0.0, rng)
    @test which in [:quick, :gradual]
end


@safetestset samplerbuilder_likelihood_cnt = "SamplerBuilder likelihood_cnt enables path_likelihood" begin
    using CompetingClocks: SamplerBuilder

    # likelihood_cnt > 1 should automatically enable path_likelihood
    builder = SamplerBuilder(Int64, Float64; likelihood_cnt=3)
    @test builder.path_likelihood == true
    @test builder.likelihood_cnt == 3
end
