using SafeTestsets

@safetestset with_common_random_smoke = "Common Random smoke" begin
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    rng = Xoshiro(234243234)
    (K, T) = (Int64, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, common_random=true)
    sampler = SamplingContext(builder, rng)
    watched = Set(1:5)
    for startup in watched
        enable!(sampler, startup, Exponential(), 0.0, 0.0)
    end
    @test length(enabled(sampler)) == 5
    time_now = 0.0
    for out in 1:5
        (when, which) = next(sampler, time_now)
        @test which ∈ watched
        fire!(sampler, which, when)
        time_now = when
        pop!(watched, which)
    end
    @test length(enabled(sampler)) == 0
end
