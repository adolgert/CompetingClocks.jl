using SafeTestsets

@safetestset with_common_random_smoke = "Common Random smoke" begin
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    rng = Xoshiro(234243234)
    (K, T) = (Int64, Float64)
    builder = SamplerBuilder(K, T; method=FirstToFireMethod(), common_random=true)
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
    builder = SamplerBuilder(K, T; method=FirstToFireMethod(), common_random=true)
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

    freeze_crn!(sampler)
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
    builder = SamplerBuilder(K, T; method=FirstToFireMethod(), common_random=true)
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

    freeze_crn!(sampler)
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
    freeze_crn!(sampler)
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


@safetestset crn_clone = "CommonRandom clone" begin
    using CompetingClocks: CommonRandom, clone, with_common_rng
    using Random: Xoshiro

    rng = Xoshiro(123456)

    # Create and populate a CommonRandom
    cr = CommonRandom{Int,Xoshiro}()

    # Record some values
    with_common_rng(cr, 1, rng) do r
        rand(r)  # Use the rng to trigger recording
    end
    with_common_rng(cr, 2, rng) do r
        rand(r)
    end

    @test haskey(cr.record, 1)
    @test haskey(cr.record, 2)

    # Clone should create an empty CommonRandom with same type parameters
    cloned = clone(cr)
    @test cloned isa CommonRandom{Int,Xoshiro}
    @test isempty(cloned.record)
    @test isempty(cloned.sample_index)
    @test isempty(cloned.miss)
    @test cloned.mode == :record
end


@safetestset crn_reset_crn = "CommonRandom reset_crn!" begin
    using CompetingClocks: CommonRandom, reset_crn!, with_common_rng
    using Random: Xoshiro

    rng = Xoshiro(789012)

    cr = CommonRandom{Int,Xoshiro}()

    # Record some values
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end
    with_common_rng(cr, 2, rng) do r
        rand(r)
    end
    with_common_rng(cr, 3, rng) do r
        rand(r)
    end

    @test length(cr.record) == 3
    @test length(cr.sample_index) == 3
    @test length(cr.miss) == 3

    # reset_crn! should clear everything
    reset_crn!(cr)
    @test isempty(cr.record)
    @test isempty(cr.sample_index)
    @test isempty(cr.miss)

    # Should be able to use it again
    with_common_rng(cr, 10, rng) do r
        rand(r)
    end
    @test haskey(cr.record, 10)
end


@safetestset crn_misscount_and_misses = "CommonRandom misscount and misses" begin
    using CompetingClocks: CommonRandom, misscount, misses, freeze_crn!, with_common_rng
    using Random: Xoshiro

    rng = Xoshiro(345678)

    cr = CommonRandom{Int,Xoshiro}()

    # Record some values for clocks 1 and 2
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end
    with_common_rng(cr, 2, rng) do r
        rand(r)
    end

    # In record mode, every call is a "miss" (no previous value)
    @test misscount(cr) == 2
    miss_pairs = collect(misses(cr))
    @test length(miss_pairs) == 2

    # Switch to replay mode
    freeze_crn!(cr)

    # Replay existing clocks (should not add to miss count since they exist)
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end
    with_common_rng(cr, 2, rng) do r
        rand(r)
    end

    # Miss count should still be 0 after freeze (freeze resets miss)
    @test misscount(cr) == 0

    # Now try to access a clock that wasn't recorded
    with_common_rng(cr, 99, rng) do r
        rand(r)
    end

    # This should add a miss
    @test misscount(cr) == 1
    miss_pairs = collect(misses(cr))
    @test any(p -> p.first == 99, miss_pairs)
end


@safetestset crn_push_to_existing_samples = "CommonRandom push to existing samples" begin
    using CompetingClocks: CommonRandom, with_common_rng, freeze_crn!, misscount
    using Random: Xoshiro

    rng = Xoshiro(901234)

    cr = CommonRandom{Int,Xoshiro}()

    # Record first sample for clock 1
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end

    @test haskey(cr.record, 1)
    @test length(cr.record[1]) == 1

    # Record second sample for same clock 1 (this exercises line 114: push to existing)
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end

    # Should now have 2 samples for clock 1
    @test length(cr.record[1]) == 2

    # Record third sample
    with_common_rng(cr, 1, rng) do r
        rand(r)
    end

    @test length(cr.record[1]) == 3

    # Now freeze and replay - all three should replay correctly
    freeze_crn!(cr)

    replay_rng = Xoshiro(999999)  # Different seed

    results = Float64[]
    for _ in 1:3
        with_common_rng(cr, 1, replay_rng) do r
            push!(results, rand(r))
        end
    end

    # The replayed values should be consistent (using saved RNG states)
    @test length(results) == 3
    @test misscount(cr) == 0  # All replays found saved values
end
