using Test
using Random
using Random: Xoshiro
using Distributions
using CompetingClocks

@testset "DelayedState: duration distribution, not pre-sampled times" begin
    rng = MersenneTwister(1234)

    # Builder with delayed support enabled
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed = true,
                             method = FirstToFireMethod())

    ctx = SamplingContext(builder, rng)

    @test Base.keytype(ctx) === Symbol
    @test timetype(ctx) === Float64
    @test ctx.delayed isa DelayedState{Symbol,Float64}

    # Enable a delayed reaction: immediate initiation, fixed duration 1.0
    delayed = Dirac(0.0) => Dirac(1.0)
    enable!(ctx, :recover, delayed)

    # The delayed state should store the *distribution* for the duration
    @test haskey(ctx.delayed.durations, :recover)
    @test ctx.delayed.durations[:recover] isa UnivariateDistribution
    @test ctx.delayed.durations[:recover] == delayed.duration

    # First event should be initiation of :recover
    when1, which1, phase1 = next_delayed(ctx)
    @test which1 === :recover
    @test phase1 === :initiate

    # Fire initiation; this should schedule completion using the stored duration distribution
    fire!(ctx, which1, phase1, when1)

    when2, which2, phase2 = next_delayed(ctx)
    @test which2 === :recover
    @test phase2 === :complete
    @test isapprox(when2 - when1, 1.0; atol=1e-8)

    # Fire completion; this should clear delayed state for :recover
    fire!(ctx, which2, phase2, when2)
    @test !haskey(ctx.delayed.durations, :recover)
end


@testset "disable! cancels all phases and clears delayed state" begin
    rng = MersenneTwister(5678)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed = true,
                             method = FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    enable!(ctx, :recover, Dirac(0.0) => Dirac(0.5))
    @test haskey(ctx.delayed.durations, :recover)
    @test isenabled(ctx, :recover)

    disable!(ctx, :recover)

    @test !isenabled(ctx, :recover)
    @test !haskey(ctx.delayed.durations, :recover)
    @test length(ctx) == 0
end


@testset "Delayed context still handles regular clocks" begin
    rng = MersenneTwister(9012)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed = true,
                             method = FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    enable!(ctx, :regular, Exponential(1.0))

    when, which, phase = next_delayed(ctx)
    @test which === :regular
    @test phase === :regular

    fire!(ctx, which, phase, when)
    @test !isenabled(ctx, :regular)
end


@testset "next_delayed errors on non-delayed context" begin
    rng = MersenneTwister(42)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed = false,
                             method = FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    @test_throws ErrorException next_delayed(ctx)

    enable!(ctx, :foo, Exponential(1.0))
    when, which = next(ctx)
    @test which === :foo
end


@testset "Re-enable during completion phase cleans up old events" begin
    rng = MersenneTwister(9012)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    # Enable first delayed reaction
    enable!(ctx, :recover, Dirac(0.0) => Dirac(5.0))

    # Fire initiation
    when1, which1, phase1 = next_delayed(ctx)
    @test phase1 === :initiate
    fire!(ctx, which1, phase1, when1)

    # Now :complete is scheduled at t=5.0
    # Re-enable with different duration
    enable!(ctx, :recover, Dirac(0.0) => Dirac(2.0))

    # Fire new initiation
    when2, which2, phase2 = next_delayed(ctx)
    @test phase2 === :initiate
    fire!(ctx, which2, phase2, when2)

    # Next should be completion at t≈2.0, not t=5.0
    when3, which3, phase3 = next_delayed(ctx)
    @test phase3 === :complete
    @test isapprox(when3, when2 + 2.0; atol=1e-8)
end


@testset "Multiple concurrent delayed clocks" begin
    rng = MersenneTwister(3456)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    enable!(ctx, :recover_a, Dirac(0.0) => Dirac(3.0))
    enable!(ctx, :recover_b, Dirac(0.0) => Dirac(1.0))
    enable!(ctx, :recover_c, Dirac(0.0) => Dirac(2.0))

    # All three should initiate at t=0
    events = Symbol[]
    for _ in 1:3
        when, which, phase = next_delayed(ctx)
        @test phase === :initiate
        push!(events, which)
        fire!(ctx, which, phase, when)
    end
    @test Set(events) == Set([:recover_a, :recover_b, :recover_c])

    # Completions should come in order: b (t=1), c (t=2), a (t=3)
    when_b, which_b, _ = next_delayed(ctx)
    @test which_b === :recover_b
    fire!(ctx, which_b, :complete, when_b)

    when_c, which_c, _ = next_delayed(ctx)
    @test which_c === :recover_c
    fire!(ctx, which_c, :complete, when_c)

    when_a, which_a, _ = next_delayed(ctx)
    @test which_a === :recover_a
    fire!(ctx, which_a, :complete, when_a)
end


@testset "clone, reset with delayed state" begin
    rng1 = MersenneTwister(7890)
    rng2 = MersenneTwister(1111)

    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod())
    ctx1 = SamplingContext(builder, rng1)

    enable!(ctx1, :recover, Dirac(0.0) => Dirac(2.0))
    fire!(ctx1, :recover, :initiate, 0.0)  # Completion now scheduled

    @test haskey(ctx1.delayed.durations, :recover)
    @test isenabled(ctx1, :recover)

    # Clone creates a fresh sampler (empty), but copies delayed state
    # Note: clone() is intended to create a "fresh" context, not copy enabled clocks
    ctx2 = clone(ctx1, rng2)
    @test haskey(ctx2.delayed.durations, :recover)
    # ctx2.sampler is fresh/empty, so the clock isn't enabled in the sampler
    # This is expected behavior - clone creates a fresh context

    # Reset
    reset!(ctx1)
    @test !haskey(ctx1.delayed.durations, :recover)
    @test !isenabled(ctx1, :recover)

    # ctx2 delayed state should be unaffected by ctx1 reset
    @test haskey(ctx2.delayed.durations, :recover)
end


@testset "CRN + Delayed: exact replay of delayed reactions" begin
    rng = Xoshiro(234243234)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod(),
                             common_random=true)
    ctx = SamplingContext(builder, rng)

    # Enable delayed reactions with stochastic durations
    enable!(ctx, :recover_a, Exponential(0.1) => Exponential(2.0))
    enable!(ctx, :recover_b, Exponential(0.2) => Exponential(1.5))

    # Record the trajectory
    trace = Tuple{Float64,Symbol,Symbol}[]
    for _ in 1:4  # 2 initiations + 2 completions
        when, which, phase = next_delayed(ctx)
        push!(trace, (when, which, phase))
        fire!(ctx, which, phase, when)
    end

    # Freeze CRN and replay
    freeze_crn!(ctx)
    enable!(ctx, :recover_a, Exponential(0.1) => Exponential(2.0))
    enable!(ctx, :recover_b, Exponential(0.2) => Exponential(1.5))

    # Replay should produce identical trajectory
    total_diff = 0.0
    for i in 1:4
        when, which, phase = next_delayed(ctx)
        @test which === trace[i][2]
        @test phase === trace[i][3]
        total_diff += abs(when - trace[i][1])
        fire!(ctx, which, phase, when)
    end

    @test total_diff < 1e-10
end


@testset "CRN + Delayed: replay with additional clocks" begin
    rng = Xoshiro(987654321)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod(),
                             common_random=true)
    ctx = SamplingContext(builder, rng)

    # First run: single delayed reaction
    enable!(ctx, :recover, Exponential(0.1) => Exponential(1.0))

    trace = Tuple{Float64,Symbol,Symbol}[]
    for _ in 1:2  # initiation + completion
        when, which, phase = next_delayed(ctx)
        push!(trace, (when, which, phase))
        fire!(ctx, which, phase, when)
    end

    # Freeze and replay with additional clock
    freeze_crn!(ctx)
    enable!(ctx, :recover, Exponential(0.1) => Exponential(1.0))
    enable!(ctx, :extra, Exponential(0.5))  # Regular clock, not delayed

    # The original delayed clock should replay exactly
    events_seen = 0
    total_diff = 0.0
    for _ in 1:3  # 2 delayed phases + 1 regular
        when, which, phase = next_delayed(ctx)
        if which === :recover
            events_seen += 1
            @test phase === trace[events_seen][3]
            total_diff += abs(when - trace[events_seen][1])
        else
            @test which === :extra
            @test phase === :regular
        end
        fire!(ctx, which, phase, when)
    end

    @test events_seen == 2
    @test total_diff < 1e-10
end


@testset "Likelihood + Delayed: steploglikelihood with delayed reactions" begin
    rng = Xoshiro(111222333)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod(),
                             path_likelihood=true)
    ctx = SamplingContext(builder, rng)

    # Enable a delayed reaction with known distributions
    enable!(ctx, :recover, Dirac(0.0) => Dirac(1.0))

    # Fire initiation at t=0
    when1, which1, phase1 = next_delayed(ctx)
    @test phase1 === :initiate
    # For delayed contexts, steploglikelihood expects internal key (clock, phase)
    ll1 = steploglikelihood(ctx, when1, (which1, phase1))
    @test isfinite(ll1)
    fire!(ctx, which1, phase1, when1)

    # Fire completion at t=1
    when2, which2, phase2 = next_delayed(ctx)
    @test phase2 === :complete
    ll2 = steploglikelihood(ctx, when2, (which2, phase2))
    @test isfinite(ll2)
    fire!(ctx, which2, phase2, when2)

    # Path likelihood should be finite
    pl = pathloglikelihood(ctx, when2)
    @test isfinite(pl)
end


@testset "Likelihood + Delayed: mixed regular and delayed clocks" begin
    rng = Xoshiro(444555666)
    builder = SamplerBuilder(Symbol, Float64;
                             support_delayed=true,
                             method=FirstToFireMethod(),
                             path_likelihood=true)
    ctx = SamplingContext(builder, rng)

    # Enable both regular and delayed clocks
    enable!(ctx, :infect, Exponential(1.0))
    enable!(ctx, :recover, Dirac(0.0) => Exponential(2.0))

    likelihoods = Float64[]
    for _ in 1:3
        when, which, phase = next_delayed(ctx)
        # For delayed contexts, steploglikelihood expects internal key (clock, phase)
        ll = steploglikelihood(ctx, when, (which, phase))
        push!(likelihoods, ll)
        fire!(ctx, which, phase, when)
    end

    # All likelihoods should be finite
    @test all(isfinite, likelihoods)

    # Path likelihood should be finite
    pl = pathloglikelihood(ctx, time(ctx))
    @test isfinite(pl)
end
