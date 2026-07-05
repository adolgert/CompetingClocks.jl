# test_gauntlet_smoke.jl - Fast, deterministic smoke test for the gauntlet
# runner. This is the ONLY gauntlet statistical-machinery test wired into the
# default CI suite (runtests.jl). It uses a tiny N and a fixed seed and asserts
# ONLY that the machinery runs end-to-end and returns well-formed verdicts. It
# deliberately does NOT assert statistical significance (that is the job of the
# real run in run_gauntlet.jl, kept out of default CI).
#
# It runs inside a @safetestset (its own module) so that its includes of the
# gauntlet source files do not collide with test_travel.jl, which includes the
# same files into Main earlier in the default suite.

using SafeTestsets

@safetestset "Gauntlet runner smoke" begin
    using Test
    using CompetingClocks
    using CompetingClocks: FirstReaction, SSA, TrackWatcher, NextReactionMethod, FirstReactionMethod
    using Random
    using Distributions
    using Graphs
    using Base.Threads
    using HypothesisTests

    include(joinpath(@__DIR__, "travel.jl"))
    using .TravelModel
    include(joinpath(@__DIR__, "generate_data.jl"))
    include(joinpath(@__DIR__, "mark_calibration.jl"))
    include(joinpath(@__DIR__, "doob_meyer.jl"))
    include(joinpath(@__DIR__, "ad_diagnostics.jl"))
    include(joinpath(@__DIR__, "anderson_darling.jl"))
    include(joinpath(@__DIR__, "running_score.jl"))
    include(joinpath(@__DIR__, "experiments.jl"))
    include(joinpath(@__DIR__, "runner.jl"))

    @testset "pvalue_verdict bands" begin
        @test pvalue_verdict(0.01) == :likely_bug
        @test pvalue_verdict(0.049) == :likely_bug
        @test pvalue_verdict(0.07) == :rerun
        @test pvalue_verdict(0.05) == :rerun
        @test pvalue_verdict(0.10) == :rerun
        @test pvalue_verdict(0.5) == :likely_correct
    end

    @testset "simple_condition" begin
        c = simple_condition()
        @test c isa TravelConfig
        @test c.memory == TravelMemory.forget
        @test c.dist == TravelRateDist.exponential
        @test c.delay == TravelRateDelay.none
    end

    @testset "run_gauntlet returns verdicts (NextReactionMethod)" begin
        # Tiny, fast, fixed-seed run. Machinery-only assertions.
        verdicts = run_gauntlet(
            NextReactionMethod();
            n_replications = 40,
            seed = 20260704,
            state_cnt = 4,
            history_steps = 3,
            verbose = false,
        )
        @test verdicts isa Vector
        @test length(verdicts) == 3
        tests = Set(v.test for v in verdicts)
        @test :doob_meyer_stepcumulant in tests
        @test :two_sample_ad_vs_firstreaction in tests
        for v in verdicts
            @test v.sampler == "NextReactionMethod"
            @test v.pvalue isa Float64
            @test 0.0 <= v.pvalue <= 1.0
            @test v.verdict in (:likely_bug, :rerun, :likely_correct)
            @test v.n > 0
        end
    end

    @testset "run_gauntlet MultiSampler composition (machinery only)" begin
        # Mixed composition: even-destination clocks -> CombinedNextReaction,
        # odd -> FirstToFire, on the complete graph so both groups get clocks.
        verdicts = run_gauntlet(
            multisampler_mixed_spec();
            n_replications = 40,
            seed = 777,
            state_cnt = 4,
            history_steps = 3,
            config = multisampler_condition(),
            verbose = false,
        )
        @test length(verdicts) == 3
        @test all(v -> v.sampler == "MultiSampler(even=>CNR,odd=>FirstToFire)", verdicts)
        @test all(v -> 0.0 <= v.pvalue <= 1.0, verdicts)
        @test all(v -> v.verdict in (:likely_bug, :rerun, :likely_correct), verdicts)
        # The chooser must actually route clocks to both groups.
        ms = multisampler_mixed_spec()(Int, Float64)
        rng = Xoshiro(1)
        enable!(ms, 1, Exponential(1.0), 0.0, 0.0, rng)
        enable!(ms, 2, Exponential(1.0), 0.0, 0.0, rng)
        @test ms.chosen[1] == :odd
        @test ms.chosen[2] == :even
        @test ms.propagator[:even] isa CompetingClocks.CombinedNextReaction
        @test ms.propagator[:odd] isa CompetingClocks.FirstToFire
    end

    @testset "run_gauntlet works for FirstReactionMethod (self-consistency)" begin
        verdicts = run_gauntlet(
            FirstReactionMethod();
            n_replications = 40,
            seed = 12345,
            state_cnt = 4,
            history_steps = 3,
            verbose = false,
        )
        @test length(verdicts) == 3
        @test all(v -> v.sampler == "FirstReactionMethod", verdicts)
    end
end
