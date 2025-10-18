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


@safetestset crn_replay_exact = "Common random numbers exact replay" begin
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    rng = Xoshiro(234243234)
    (K, T) = (Int64, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, common_random=true)
    sampler = SamplingContext(builder, rng)

    time_now = zero(T)
    trace = Dict{Int,Float64}()
    watched = Set(1:5)
    for startup in watched
        enable!(sampler, startup, Exponential(), 0.0, 0.0)
    end
    @test length(enabled(sampler)) == 5
    for out in 1:5
        (when, which) = next(sampler, time_now)
        trace[which] = when
        @test which ∈ watched
        fire!(sampler, which, when)
        time_now = when
        pop!(watched, which)
    end

    freeze!(sampler)
    for startup in Set(1:5)
        enable!(sampler, startup, Exponential(), 0.0, 0.0)
    end
    @test length(enabled(sampler)) == 5
    time_now = zero(T)
    total_diff::T = 0
    for out in 1:5
        (when, which) = next(sampler, time_now)
        total_diff += abs(when - trace[which])
        fire!(sampler, which, when)
        time_now = when
    end

    @test total_diff < 1e-10
end


@safetestset crn_replay_and_more = "Common random numbers replay plus some" begin
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    rng = Xoshiro(234243234)
    (K, T) = (Int64, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, common_random=true)
    sampler = SamplingContext(builder, rng)

    time_now = zero(T)
    trace = Dict{Int,Float64}()
    watched = Set(1:5)
    for startup in watched
        enable!(sampler, startup, Exponential(), 0.0, 0.0)
    end
    @test length(enabled(sampler)) == 5
    for out in 1:5
        (when, which) = next(sampler, time_now)
        trace[which] = when
        @test which ∈ watched
        fire!(sampler, which, when)
        time_now = when
        pop!(watched, which)
    end

    freeze!(sampler)
    for startup in Set(1:10)
        enable!(sampler, startup, Exponential(), 0.0, 0.0)
    end
    @test length(enabled(sampler)) == 10
    time_now = zero(T)
    total_diff::T = 0
    for out in 1:10
        (when, which) = next(sampler, time_now)
        if which <= 5
            total_diff += abs(when - trace[which])
        else
            @test isfinite(when)
        end
        fire!(sampler, which, when)
        time_now = when
    end

    @test total_diff < 1e-10
end
