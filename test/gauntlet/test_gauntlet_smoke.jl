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
    include(joinpath(@__DIR__, "config.jl"))
    include(joinpath(@__DIR__, "matrix.jl"))
    include(joinpath(@__DIR__, "reporting.jl"))

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
        enable!(ms, 1, Exponential(1.0), 0.0, 0.0)
        enable!(ms, 2, Exponential(1.0), 0.0, 0.0)
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

    # --- Matrix machinery (config.jl / matrix.jl / reporting.jl) -----------
    # Tiny N, tiny n, fixed seeds. Machinery-only: asserts well-formed TestRuns
    # and report/CSV generation, NOT statistical significance.
    @testset "config design + filtering" begin
        design = generate_full_design(; state_cnt=4, history_steps=3, doob_steps=20)
        @test length(design) == 12                      # 3 families x 2 memory x 2 graph
        @test all(c -> c isa TestConfiguration, design)
        expo = exponential_only_configs(design)
        @test length(expo) == 4                         # only :exponential family
        @test all(c -> c.family == :exponential, expo)
        @test filter_for_sampler(design, :general) == design
        @test filter_for_sampler(design, :exponential) == expo
        # round-trip key parsing
        c = design[1]
        @test parse_config_key(config_key(c); state_cnt=4, history_steps=3, doob_steps=20).family == c.family
        # registry capabilities
        @test sampler_capability("RejectionMethod(RSSA)") == :exponential
        @test sampler_capability("NextReactionMethod") == :general
    end

    @testset "run_single_cell returns TestRuns" begin
        cfg = TestConfiguration(:weibull, :forget, :cycle; state_cnt=4, history_steps=3, doob_steps=20)
        runs = run_single_cell("NextReactionMethod", NextReactionMethod(), cfg, 30, 20260705)
        @test length(runs) == 3
        @test all(r -> r isa TestRun, runs)
        @test all(r -> 0.0 <= r.pvalue <= 1.0, runs)
        @test all(r -> r.error_message === nothing, runs)
        # exponential-only sampler on a Weibull condition errors gracefully
        bad = run_single_cell("RejectionMethod(RSSA)", RejectionMethod(), cfg, 30, 1)
        @test all(r -> r.verdict == :error && isnan(r.pvalue), bad)
    end

    @testset "run_matrix + BH + reporting" begin
        # A minimal 2-sampler x exponential-only design.
        design = generate_full_design(; families=(:exponential,), memories=(:forget,),
            graphs=(:cycle, :complete), state_cnt=4, history_steps=3, doob_steps=20)
        samplers = Tuple{String,Any,Symbol}[
            ("FirstReactionMethod", FirstReactionMethod(), :general),
            ("DirectMethod(:remove,:tree)", DirectMethod(:remove, :tree), :exponential),
        ]
        runs = run_matrix(; n=30, base_seed=20260705, samplers=samplers, design=design, verbose=false)
        @test length(runs) == 2 * 2 * 3               # 2 samplers x 2 conditions x 3 tests
        @test all(r -> r isa TestRun, runs)

        adj, flagged = apply_bh(runs; alpha=0.05)
        @test length(adj) == length(runs)
        @test all(a -> isnan(a) || (0.0 <= a <= 1.0), adj)
        @test flagged isa Vector{Int}

        # BH is monotone-ish: adjusted >= raw for the smallest p-value.
        raw = [r.pvalue for r in runs]
        i = argmin(raw)
        @test adj[i] >= raw[i] - 1e-9

        # grouping helpers
        g = group_by_parameter(runs, :sampler)
        @test haskey(g, Symbol("FirstReactionMethod"))
        prob = identify_problematic_configs(runs; threshold=1.0)  # everything qualifies
        @test !isempty(prob)

        # report + CSV write to a temp dir
        d = mktempdir()
        rp = generate_markdown_report(runs, joinpath(d, "r.md"); n=30, base_seed=20260705)
        cp = write_matrix_csv(runs, joinpath(d, "r.csv"))
        @test isfile(rp) && filesize(rp) > 0
        @test isfile(cp) && filesize(cp) > 0
        @test occursin("BH FDR", read(rp, String))
    end
end
