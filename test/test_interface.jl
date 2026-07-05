using SafeTestsets


@safetestset interface_error_fallbacks = "SSA interface error fallbacks" begin
    using CompetingClocks: SSA, enable!, disable!, fire!, next, reset!, copy_clocks!,
        clone, jitter!, enabled
    using Random: Xoshiro
    using Distributions: Exponential

    # A minimal sampler that doesn't implement any methods
    struct UnimplementedSampler <: SSA{Int,Float64} end

    rng = Xoshiro(123456)
    sampler = UnimplementedSampler()

    @test_throws ErrorException enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    @test_throws ErrorException disable!(sampler, 1, 0.0)
    # fire! has a documented interface fallback that forwards to disable!,
    # so on an unimplemented sampler it surfaces disable!'s error.
    @test_throws ErrorException fire!(sampler, 1, 0.0)
    @test_throws ErrorException next(sampler, 0.0, rng)
    @test_throws ErrorException reset!(sampler)
    @test_throws ErrorException copy_clocks!(sampler, sampler)
    @test_throws ErrorException clone(sampler)
    @test_throws ErrorException jitter!(sampler, 0.0, rng)
    @test_throws ErrorException sampler[1]
    @test_throws ErrorException keys(sampler)
    @test_throws ErrorException length(sampler)
    @test_throws ErrorException haskey(sampler, 1)
    @test_throws ErrorException enabled(sampler)
end
