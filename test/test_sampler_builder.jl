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
