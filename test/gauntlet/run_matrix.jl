# run_matrix.jl - Driver for the full gauntlet statistical test matrix.
#
# NOT part of the default CI suite. Runs the whole sampler x condition matrix at
# a real sample size, applies Benjamini-Hochberg FDR correction across the
# matrix, runs the null-control seed sweep on any flagged cell (and on the
# milestone-1 cycle cell for documentation), and writes:
#   results/gauntlet/matrix_report.md   (markdown report; headline = BH verdict)
#   results/gauntlet/matrix_results.csv (every raw result + BH-adjusted p)
#
# Usage:
#   # full matrix (default n=1000, base_seed=20260705):
#   julia --project=test test/gauntlet/run_matrix.jl
#   julia --project=test test/gauntlet/run_matrix.jl 1000 20260705
#
#   # reproduce ONE cell (all three tests) at its own seed:
#   julia --project=test test/gauntlet/run_matrix.jl cell "NextReactionMethod" "exponential/forget/cycle" 1000 20260705

using CompetingClocks
using CompetingClocks: FirstReaction, SSA, TrackWatcher, NextReactionMethod,
    FirstReactionMethod, FirstToFireMethod, DirectMethod, RejectionMethod,
    PartialPropensityMethod
using Random
using Distributions
using Graphs
using Base.Threads
using HypothesisTests

include("travel.jl")
using .TravelModel
include("generate_data.jl")
include("mark_calibration.jl")
include("doob_meyer.jl")
include("ad_diagnostics.jl")
include("anderson_darling.jl")
include("running_score.jl")
include("runner.jl")
include("config.jl")
include("matrix.jl")
include("reporting.jl")


function run_full_matrix(n::Int, base_seed::Int)
    outdir = joinpath(@__DIR__, "..", "..", "results", "gauntlet")
    mkpath(outdir)
    report_path = joinpath(outdir, "matrix_report.md")
    csv_path = joinpath(outdir, "matrix_results.csv")

    t0 = time()
    runs = run_matrix(; n=n, base_seed=base_seed, verbose=true)
    println("\nMatrix compute time: $(round(time() - t0, digits=1)) s")

    adjusted, flagged = apply_bh(runs; alpha=0.05)
    print_summary_table(runs)

    # Null-control discipline: investigate every BH-flagged cell, the milestone-1
    # cell (NextReactionMethod @ exponential/forget/cycle, which showed a
    # single-seed raw two-sample flag p=0.0183 at seed 20260704), AND the
    # worst cell of any sampler whose raw sub-0.05 rate exceeds the trusted
    # FirstReactionMethod reference's own rate (so every "above reference"
    # sampler in the report gets a seed sweep).
    investigate = Set{Tuple{String,String}}()
    for i in flagged
        push!(investigate, (runs[i].sampler_name, runs[i].config_key))
    end
    push!(investigate, ("NextReactionMethod", "exponential/forget/cycle"))

    by_sampler = group_by_parameter(runs, :sampler)
    fr_rate = by_sampler[Symbol("FirstReactionMethod")].frac_below_05
    for (skey, _spec, _cap) in matrix_samplers()
        s = by_sampler[Symbol(skey)]
        if s.frac_below_05 > fr_rate + 1e-9
            # pick this sampler's lowest-p cell
            cells = [r for r in runs if r.sampler_name == skey && !isnan(r.pvalue)]
            worst_cell = cells[argmin([r.pvalue for r in cells])]
            push!(investigate, (skey, worst_cell.config_key))
        end
    end

    investigations = NamedTuple[]
    sweep_seeds = base_seed .+ (0:9) .* 101   # 10 well-separated seeds
    for (skey, ckey) in sort(collect(investigate))
        println("\nNull-control investigation: $skey @ $ckey")
        cfg = parse_config_key(ckey; doob_steps=n)
        inv = investigate_flag(skey, cfg; seeds=sweep_seeds, n=n)
        push!(investigations, inv)
        println("  SUT  frac p<0.05: DM=$(round(inv.sut_frac.doob,digits=2)) ",
            "AD=$(round(inv.sut_frac.two_sample,digits=2)) MK=$(round(inv.sut_frac.mark,digits=2))")
        println("  null frac p<0.05: DM=$(round(inv.null_frac.doob,digits=2)) ",
            "AD=$(round(inv.null_frac.two_sample,digits=2)) MK=$(round(inv.null_frac.mark,digits=2))")
    end

    write_matrix_csv(runs, csv_path; adjusted=adjusted)
    generate_markdown_report(runs, report_path; n=n, base_seed=base_seed,
        alpha=0.05, investigations=investigations)
    println("\nWrote:\n  $report_path\n  $csv_path")
    return runs
end


function run_one_cell(sampler_key::String, config_key_str::String, n::Int, seed::Int)
    spec = sampler_spec_for(sampler_key)
    cfg = parse_config_key(config_key_str; doob_steps=n)
    println("Reproducing cell: $sampler_key @ $config_key_str  (n=$n, seed=$seed)")
    runs = run_single_cell(sampler_key, spec, cfg, n, seed)
    println(rpad("test", 34), lpad("p-value", 10), "  verdict")
    println("-"^60)
    for r in runs
        println(rpad(string(r.test_type), 34), lpad(round(r.pvalue, digits=4), 10),
            "  ", r.verdict, r.error_message === nothing ? "" : "  ($(r.error_message))")
    end
    return runs
end


function main()
    if length(ARGS) >= 1 && ARGS[1] == "cell"
        length(ARGS) >= 3 || error("usage: run_matrix.jl cell <sampler_key> <config_key> [n] [seed]")
        sampler_key = ARGS[2]
        config_key_str = ARGS[3]
        n = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 1000
        seed = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 20260705
        run_one_cell(sampler_key, config_key_str, n, seed)
    else
        n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
        base_seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260705
        run_full_matrix(n, base_seed)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
