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
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 5
    for out in 1:5
        (when, which) = next(sampler)
        @test which ∈ watched
        fire!(sampler, which, when)
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

    trace = Dict{Int,Float64}()
    watched = Set(1:5)
    for startup in watched
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 5
    for out in 1:5
        (when, which) = next(sampler)
        trace[which] = when
        @test which ∈ watched
        fire!(sampler, which, when)
        pop!(watched, which)
    end

    freeze!(sampler)
    for startup in Set(1:5)
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 5
    total_diff::T = 0
    for out in 1:5
        (when, which) = next(sampler)
        total_diff += abs(when - trace[which])
        fire!(sampler, which, when)
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

    trace = Dict{Int,Float64}()
    watched = Set(1:5)
    for startup in watched
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 5
    for out in 1:5
        (when, which) = next(sampler)
        trace[which] = when
        @test which ∈ watched
        fire!(sampler, which, when)
        pop!(watched, which)
    end

    freeze!(sampler)
    for startup in Set(1:10)
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 10
    total_diff::T = 0
    extras = Dict{Int,T}()
    for out in 1:10
        (when, which) = next(sampler)
        if which <= 5
            total_diff += abs(when - trace[which])
        else
            @test isfinite(when)
            extras[which] = when
        end
        fire!(sampler, which, when)
    end
    @test total_diff < 1e-10

    # Try again. Check that the last 5 remain random.
    freeze!(sampler)
    for startup in Set(1:10)
        enable!(sampler, startup, Exponential())
    end
    @test length(enabled(sampler)) == 10
    total_diff = zero(T)
    extras_diff = zero(T)
    for out in 1:10
        (when, which) = next(sampler)
        if which <= 5
            total_diff += abs(when - trace[which])
        else
            extras_diff += abs(when - extras[which])
        end
        fire!(sampler, which, when)
    end

    @test total_diff < 1e-10
    @test extras_diff > 0.1
end
