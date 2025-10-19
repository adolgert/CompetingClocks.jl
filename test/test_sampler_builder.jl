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
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    @test isempty(scbuild.group)
    add_group!(scbuild, :capistrano => (k, d) -> true; sampler_spec=:firsttofire)
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_once = "SamplerBuilder single construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple, Float64; sampler_spec=:firsttofire)
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    @test !isempty(scbuild.group)
    scresult = build_sampler(scbuild)
    @test scresult isa FirstToFire
end

@safetestset samplerbuilder_construct_multi = "SamplerBuilder multi construction" begin
    using CompetingClocks
    scbuild = SamplerBuilder(Tuple, Float64)
    @test (:direct, :keep, :tree) in available_samplers(scbuild)
    add_group!(scbuild, :sparky => (x, d) -> x[1] == :recover, sampler_spec=(:nextreaction,))
    add_group!(scbuild, :forthright => (x, d) -> x[1] == :infect)
    scresult = build_sampler(scbuild)
    @test scresult isa MultiSampler
end
