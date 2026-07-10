# For this macro, need to import it into main for it to work.
import InteractiveUtils: @code_typed



@testset "SamplerContext doesn't penalize missing pieces" begin
    using CompetingClocks
    using Random
    using Distributions

    K = Int64
    T = Float64
    rng = Xoshiro(899987987)
    builder = SamplerBuilder(K, T; method=FirstToFireMethod())
    sampler = SamplingContext(builder, rng)
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity))
    end
    # Assert that when compiled, the if-then statements in the context
    # are compiled away!
    res = @code_typed enable!(sampler, 1, Exponential(0.5), 0.0)
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
    for SamplerType in [FirstToFireMethod(), DirectMethod(), FirstReactionMethod(), NextReactionMethod()]
        rng = Xoshiro(90422342)
        builder = SamplerBuilder(K, T; method=SamplerType)
        context = SamplingContext(builder, rng)
        test_enable = Set{Int64}()
        for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
            enable!(context, clock_id, Exponential(propensity))
            push!(test_enable, clock_id)
        end
        when, which = next(context)
        fire!(context, which, when)
        delete!(test_enable, which)
        when, which = next(context)
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


@safetestset context_short_constructor = "SamplingContext short constructor" begin
    using CompetingClocks: SamplingContext, enable!, next, fire!, keytype, timetype
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(123456)
    # Use the short constructor SamplingContext(K, T, rng; kwargs...)
    ctx = SamplingContext(Int64, Float64, rng)

    @test keytype(ctx) == Int64
    @test timetype(ctx) == Float64

    enable!(ctx, 1, Exponential(1.0))
    when, which = next(ctx)
    @test which == 1
    @test when > 0.0
end


@safetestset context_with_likelihood_cnt = "SamplingContext with likelihood_cnt > 1" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, next, fire!
    using CompetingClocks: sample_from_distribution!, pathloglikelihood
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(234567)
    # likelihood_cnt > 1 creates PathLikelihoods
    builder = SamplerBuilder(Int64, Float64; likelihood_cnt=2)
    ctx = SamplingContext(builder, rng)

    # Enable with vector of distributions
    enable!(ctx, 1, [Exponential(1.0), Exponential(2.0)], 0.0)
    enable!(ctx, 2, [Exponential(0.5), Exponential(1.0)], 0.0)

    # Can change which distribution to sample from
    sample_from_distribution!(ctx, 2)
    @test ctx.sample_distribution == 2

    when, which = next(ctx)
    fire!(ctx, which, when)

    # pathloglikelihood returns a vector when likelihood_cnt > 1
    ll = pathloglikelihood(ctx, 1.0)
    @test length(ll) == 2
end


@safetestset context_with_step_likelihood = "SamplingContext with step_likelihood" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, next, fire!
    using CompetingClocks: steploglikelihood, FirstToFireMethod, enabled_history, disabled_history
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(345678)
    # step_likelihood with a sampler that doesn't have built-in step likelihood
    # should create a TrackWatcher
    builder = SamplerBuilder(Int64, Float64; step_likelihood=true, method=FirstToFireMethod())
    ctx = SamplingContext(builder, rng)

    enable!(ctx, 1, Exponential(1.0))
    enable!(ctx, 2, Exponential(2.0))

    when, which = next(ctx)

    # steploglikelihood should work via the likelihood object
    ll = steploglikelihood(ctx, when, which)
    @test isfinite(ll)

    fire!(ctx, which, when)

    @test_throws "`recording=true`" enabled_history(ctx)
    @test_throws "`recording=true`" disabled_history(ctx)
end


@safetestset context_with_debug = "SamplingContext with debug/recording" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, next, fire!, disable!
    using CompetingClocks: enabled_history, disabled_history, debug_of
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    # recording=true creates DebugWatcher
    builder = SamplerBuilder(Int64, Float64; recording=true)
    ctx = SamplingContext(builder, rng)

    @test debug_of(ctx) !== nothing

    enable!(ctx, 1, Exponential(1.0))
    enable!(ctx, 2, Exponential(2.0))

    when, which = next(ctx)
    fire!(ctx, which, when)
    @test length(disabled_history(ctx)) == 1

    # Disable the other clock
    other = which == 1 ? 2 : 1
    disable!(ctx, other)

    @test length(ctx) == 0
    @test length(enabled_history(ctx)) == 2
    @test length(disabled_history(ctx)) == 2
end


@safetestset context_clone = "SamplingContext clone" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, clone
    using Random: Xoshiro
    using Distributions: Exponential

    rng1 = Xoshiro(567890)
    rng2 = Xoshiro(678901)

    builder = SamplerBuilder(Int64, Float64; path_likelihood=true)
    ctx = SamplingContext(builder, rng1)
    enable!(ctx, 1, Exponential(1.0))

    # Clone with a different RNG
    cloned = clone(ctx, rng2)

    @test cloned.rng === rng2
    @test cloned.sampler !== ctx.sampler  # Different sampler instance
    @test length(cloned) == 0  # Cloned sampler is empty
    @test length(ctx) == 1  # Original unchanged
end


@safetestset context_reset = "SamplingContext reset!" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, next, fire!, reset!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(789012)
    builder = SamplerBuilder(Int64, Float64; path_likelihood=true)
    ctx = SamplingContext(builder, rng)

    enable!(ctx, 1, Exponential(1.0))
    enable!(ctx, 2, Exponential(2.0))
    when, which = next(ctx)
    fire!(ctx, which, when)

    @test time(ctx) > 0.0

    reset!(ctx)

    @test time(ctx) == 0.0
    @test length(ctx) == 0
end


@safetestset context_vector_enable = "SamplingContext enable! with vector distribution" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, next, fire!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(890123)
    builder = SamplerBuilder(Int64, Float64; likelihood_cnt=2)
    ctx = SamplingContext(builder, rng)

    # Enable with vector of distributions (exercises lines 271-285)
    enable!(ctx, 1, [Exponential(1.0), Exponential(2.0)], 0.0)

    # Enable with vector but no relative_te (exercises lines 304-305)
    enable!(ctx, 2, [Exponential(0.5), Exponential(1.0)])

    @test length(ctx) == 2

    when, which = next(ctx)
    @test which in [1, 2]
end


@safetestset context_length_with_likelihood = "SamplingContext length/enabled with likelihood" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, enable!, enabled, isenabled
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(901234)
    # path_likelihood creates a TrajectoryWatcher which tracks enabled clocks
    builder = SamplerBuilder(Int64, Float64; path_likelihood=true)
    ctx = SamplingContext(builder, rng)

    enable!(ctx, 1, Exponential(1.0))
    enable!(ctx, 2, Exponential(2.0))

    # These should delegate to the likelihood object
    @test length(ctx) == 2
    @test 1 in enabled(ctx)
    @test 2 in enabled(ctx)
    @test isenabled(ctx, 1)
    @test isenabled(ctx, 2)
    @test !isenabled(ctx, 999)
end


@safetestset context_sample_from_distribution_errors = "SamplingContext sample_from_distribution! errors" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, sample_from_distribution!
    using Random: Xoshiro

    rng = Xoshiro(012345)

    # Without likelihood_cnt > 1, can't sample from later distribution
    builder1 = SamplerBuilder(Int64, Float64)
    ctx1 = SamplingContext(builder1, rng)
    @test_throws ErrorException sample_from_distribution!(ctx1, 2)

    # With likelihood_cnt, can't sample from out-of-range index
    builder2 = SamplerBuilder(Int64, Float64; likelihood_cnt=3)
    ctx2 = SamplingContext(builder2, rng)
    @test_throws ErrorException sample_from_distribution!(ctx2, 5)  # out of range
end


# The old common-random-numbers runtime (freeze_crn! / reset_crn! on a
# common_random=true context) is retired. Common random numbers are now a
# property of the sampler's keyed streams: build two contexts from the same seed
# and they draw identically per clock. That successor is exercised in
# test_crn_streams.jl.


@safetestset context_likelihood_capability_probe = "SamplingContext likelihood capability probes" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, DirectMethod, FirstToFireMethod
    using CompetingClocks: enable!, next, steploglikelihood, pathloglikelihood, likelihood_of
    using Random: Xoshiro
    using Distributions: Exponential

    # (1) A context built with path_likelihood=true still returns pathloglikelihood
    #     (existing behavior, via the likelihood watcher).
    rng = Xoshiro(20240704)
    builder_path = SamplerBuilder(Int64, Float64; path_likelihood=true)
    ctx_path = SamplingContext(builder_path, rng)
    enable!(ctx_path, 1, Exponential(1.0))
    next(ctx_path)
    @test pathloglikelihood(ctx_path, 1.0) isa Real

    # (2) A context whose sampler lacks the capability (FirstToFire supports
    #     neither step nor path likelihood) and has no likelihood watcher throws
    #     the informative "doesn't support" error.
    builder_ff = SamplerBuilder(Int64, Float64; method=FirstToFireMethod())
    ctx_ff = SamplingContext(builder_ff, rng)
    @test likelihood_of(ctx_ff) === nothing
    enable!(ctx_ff, 1, Exponential(1.0))
    next(ctx_ff)
    @test_throws ErrorException pathloglikelihood(ctx_ff, 1.0)
    @test_throws ErrorException steploglikelihood(ctx_ff, 1.0, 1)

    # (3) Regression: errors from within the sampler's likelihood code must
    #     propagate, NOT be swallowed into the misleading "doesn't support"
    #     message. A bare DirectCall-backed context (trait true, no watcher)
    #     querying a nonexistent clock key raises a KeyError.
    builder_dc = SamplerBuilder(Int64, Float64; method=DirectMethod())
    ctx_dc = SamplingContext(builder_dc, rng)
    @test likelihood_of(ctx_dc) === nothing
    enable!(ctx_dc, 1, Exponential(1.0))
    next(ctx_dc)
    # A valid key returns a real value (trait routes directly to the sampler).
    @test steploglikelihood(ctx_dc, 1.0, 1) isa Real
    @test pathloglikelihood(ctx_dc, 1.0) isa Real
    # A nonexistent key surfaces the real underlying error, not the misleading one.
    @test_throws KeyError steploglikelihood(ctx_dc, 1.0, 999999)
end


@safetestset context_custom_watcher = "SamplingContext one-place extensibility via watcher tuple" begin
    # Demonstrates that adding a brand-new observer watcher requires touching
    # ONLY the watcher tuple: the context fan-out (context.jl) is not modified.
    using CompetingClocks
    using CompetingClocks: SamplingContext, next, FirstToFire
    using Random: Xoshiro
    using Distributions: Exponential

    # A ~10-line custom watcher counting the three lifecycle callbacks.
    mutable struct CountingWatcher
        enables::Int
        disables::Int
        fires::Int
    end
    CountingWatcher() = CountingWatcher(0, 0, 0)
    CompetingClocks.enable!(w::CountingWatcher, clock, dist, te, when) = (w.enables += 1; nothing)
    CompetingClocks.disable!(w::CountingWatcher, clock, when) = (w.disables += 1; nothing)
    CompetingClocks.fire!(w::CountingWatcher, clock, when) = (w.fires += 1; nothing)
    CompetingClocks.reset!(w::CountingWatcher) = (w.enables = w.disables = w.fires = 0; nothing)
    CompetingClocks.clone(w::CountingWatcher) = CountingWatcher()

    rng = Xoshiro(2024)
    sampler = FirstToFire{Int64,Float64}()
    watcher = CountingWatcher()
    # The ONLY thing needed to support the new observer: put it in the tuple.
    ctx = SamplingContext{Int64,Float64,typeof(sampler),typeof(rng),
                          Tuple{CountingWatcher},Nothing}(
        sampler, rng, (watcher,), 1.0, 0.0, 0.0, 1, nothing)

    enable!(ctx, 1, Exponential(1.0))
    enable!(ctx, 2, Exponential(2.0))
    @test watcher.enables == 2

    when, which = next(ctx)
    fire!(ctx, which, when)
    @test watcher.fires == 1

    other = which == 1 ? 2 : 1
    disable!(ctx, other)
    @test watcher.disables == 1

    reset!(ctx)
    @test watcher.enables == 0
    @test watcher.disables == 0
    @test watcher.fires == 0
end
