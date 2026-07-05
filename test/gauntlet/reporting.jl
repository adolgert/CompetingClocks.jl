# reporting.jl - Phase 4 of the gauntlet: reporting over a Vector{TestRun}.
#
# Produces:
#   - a markdown report (matrix verdict table + flagged-cell detail + summaries)
#   - a CSV of every raw result (with BH-adjusted p-values)
#   - the interactive-analysis helpers `identify_problematic_configs` and
#     `group_by_parameter`.
#
# The headline verdict is the Benjamini-Hochberg FDR correction applied ACROSS
# the whole matrix (see matrix.jl `apply_bh`); raw per-cell p-values are always
# retained alongside it.
#
# Assumes config.jl and matrix.jl are in scope.


const TEST_SHORT = Dict(
    :doob_meyer_stepcumulant => "DM",
    :two_sample_ad_vs_firstreaction => "AD",
    :mark_calibration => "MK",
)

const VERDICT_GLYPH = Dict(
    :likely_correct => ".",
    :rerun => "~",
    :likely_bug => "X",
    :error => "E",
)


"""
    repro_command(run; n) -> String

The exact shell command that reproduces a single cell (all three tests) at the
cell's own seed.
"""
function repro_command(run::TestRun; n::Int=1000)
    return string("julia --project=test test/gauntlet/run_matrix.jl cell ",
        "\"", run.sampler_name, "\" \"", run.config_key, "\" ", n, " ", run.seed)
end


"""
    identify_problematic_configs(runs; threshold) -> Vector

Configurations (as (sampler, config_key) pairs) with at least one raw p-value
below `threshold`. Sorted by the minimum raw p-value.
"""
function identify_problematic_configs(runs::Vector{TestRun}; threshold::Float64=0.05)
    groups = Dict{Tuple{String,String},Vector{TestRun}}()
    for r in runs
        push!(get!(groups, (r.sampler_name, r.config_key), TestRun[]), r)
    end
    bad = [(k, minimum(x -> isnan(x.pvalue) ? Inf : x.pvalue, v))
           for (k, v) in groups
           if any(x -> !isnan(x.pvalue) && x.pvalue < threshold, v)]
    sort!(bad, by=x -> x[2])
    return bad
end


"""
    group_by_parameter(runs, param) -> Dict

Group runs by a condition parameter (`:family`, `:memory`, `:graph`, `:sampler`,
or `:test`) and report each group's below-0.05 rate. Reveals whether flags
cluster on one axis (a pattern) rather than scattering (multiple-testing noise).
"""
function group_by_parameter(runs::Vector{TestRun}, param::Symbol)
    keyof(r) =
        param == :family ? r.config.family :
        param == :memory ? r.config.memory :
        param == :graph ? r.config.graph :
        param == :sampler ? Symbol(r.sampler_name) :
        param == :test ? r.test_type :
        error("unknown param: $param")
    groups = Dict{Any,Vector{TestRun}}()
    for r in runs
        push!(get!(groups, keyof(r), TestRun[]), r)
    end
    summary = Dict{Any,NamedTuple}()
    for (k, v) in groups
        valid = filter(x -> !isnan(x.pvalue), v)
        nbelow = count(x -> x.pvalue < 0.05, valid)
        summary[k] = (; n=length(valid),
            frac_below_05=isempty(valid) ? 0.0 : nbelow / length(valid),
            min_p=isempty(valid) ? NaN : minimum(x -> x.pvalue, valid))
    end
    return summary
end


"""
    write_matrix_csv(runs, path; adjusted)

Write every raw result to CSV, including the BH-adjusted p-value.
"""
function write_matrix_csv(runs::Vector{TestRun}, path::AbstractString;
        adjusted::Union{Nothing,Vector{Float64}}=nothing)
    adj = adjusted === nothing ? benjamini_hochberg([r.pvalue for r in runs]) : adjusted
    open(path, "w") do io
        println(io, "sampler,family,memory,graph,config_key,test,seed,pvalue,bh_pvalue,verdict,n,elapsed_s,error")
        for (i, r) in enumerate(runs)
            err = r.error_message === nothing ? "" : replace(r.error_message, ","=>";", "\n"=>" ")
            println(io, join((
                "\"" * r.sampler_name * "\"",
                r.config.family, r.config.memory, r.config.graph,
                r.config_key, r.test_type, r.seed,
                isnan(r.pvalue) ? "NA" : round(r.pvalue, digits=6),
                isnan(adj[i]) ? "NA" : round(adj[i], digits=6),
                r.verdict, r.n, round(r.elapsed, digits=3),
                "\"" * err * "\"",
            ), ","))
        end
    end
    return path
end


# --- Markdown report --------------------------------------------------------

function _cell_triple(runs_by_test)
    # runs_by_test: Dict test_type => TestRun for one (sampler,config)
    parts = String[]
    for t in (:doob_meyer_stepcumulant, :two_sample_ad_vs_firstreaction, :mark_calibration)
        r = get(runs_by_test, t, nothing)
        push!(parts, r === nothing ? "-" : VERDICT_GLYPH[r.verdict])
    end
    return join(parts, "")
end


"""
    generate_markdown_report(runs, path; n, base_seed, alpha)

Write the full matrix report to `path`. Sections: executive summary with the BH
headline, the sampler x condition verdict matrix (each cell a DM/AD/MK verdict
triple), flagged-cell detail (raw p, BH-adjusted p, repro command), and
per-parameter grouping summaries.
"""
function generate_markdown_report(runs::Vector{TestRun}, path::AbstractString;
        n::Int=1000, base_seed::Integer=20260705, alpha::Float64=0.05,
        investigations::Vector=NamedTuple[])
    adjusted, flagged = apply_bh(runs; alpha=alpha)
    nerr = count(r -> r.verdict == :error, runs)
    nbug = count(r -> r.verdict == :likely_bug, runs)
    nrerun = count(r -> r.verdict == :rerun, runs)
    valid = count(r -> !isnan(r.pvalue), runs)

    # index by (sampler, config_key)
    samplers = matrix_samplers()
    sampler_order = [k for (k, _s, _c) in samplers]
    design = generate_full_design(; doob_steps=n)
    col_keys = [config_key(c) for c in design]

    by_sc = Dict{Tuple{String,String},Dict{Symbol,TestRun}}()
    for r in runs
        d = get!(by_sc, (r.sampler_name, r.config_key), Dict{Symbol,TestRun}())
        d[r.test_type] = r
    end

    open(path, "w") do io
        println(io, "# Gauntlet Statistical Test-Matrix Report")
        println(io)
        println(io, "Reproducible from `(base_seed=$(base_seed), n=$(n))`. ",
            "Fixed seeds; each cell seed is `base_seed + 10*cell_index` in registry x design order.")
        println(io)
        println(io, "## Executive Summary")
        println(io)
        println(io, "- Total results: **$(length(runs))** ($(length(runs) ÷ 3) cells x 3 tests), valid p-values: **$valid**, errors: **$nerr**")
        println(io, "- Raw per-cell bands: `likely_bug` (p<0.05): **$nbug**, `rerun` (0.05<=p<=0.10): **$nrerun**")
        println(io, "- **BH FDR headline (alpha=$(alpha) across the whole matrix): $(length(flagged)) result(s) flagged.**")
        if isempty(flagged)
            println(io)
            println(io, "> No cell survives Benjamini-Hochberg correction across the matrix. ",
                "Raw sub-0.05 p-values are consistent with the false-positive rate expected from ",
                "$(valid) simultaneous tests. No sampler shows a statistically-supported defect.")
        end
        println(io)

        # Automatic null-control interpretation: the trusted FirstReactionMethod
        # reference, run as its own SUT through the identical path, is the built-in
        # null. If no SUT exceeds its raw sub-0.05 rate, the flags are noise.
        by_sampler = group_by_parameter(runs, :sampler)
        fr_key = Symbol("FirstReactionMethod")
        if haskey(by_sampler, fr_key)
            fr_rate = by_sampler[fr_key].frac_below_05
            worst = maximum(s.frac_below_05 for s in values(by_sampler))
            worst_names = sort([string(k) for (k, s) in by_sampler if s.frac_below_05 >= worst - 1e-9])
            println(io, "**Null-control discipline.** The trusted `FirstReactionMethod` reference is ",
                "run as its own SUT (FirstReaction-vs-FirstReaction two-sample, plus its own ",
                "Doob-Meyer/Mark), a built-in null. Its raw sub-0.05 rate is ",
                "**$(round(fr_rate, digits=3))**. The maximum sub-0.05 rate over all samplers is ",
                "**$(round(worst, digits=3))** (", join(worst_names, ", "), "). ",
                fr_rate >= worst - 1e-9 ?
                    "No sampler under test exceeds the trusted reference's own false-positive rate, so the raw flags are attributable to multiple-testing noise, not sampler defects." :
                    "Samplers above the reference rate are examined in the null-control investigations below.")
            println(io)
        end
        println(io, "Verdict glyphs: `.`=likely_correct (p>0.10), `~`=rerun (0.05-0.10), ",
            "`X`=likely_bug (p<0.05, raw), `E`=error, `-`=not run (condition excluded for this sampler). ",
            "Each cell is a **DM/AD/MK** triple: Doob-Meyer / two-sample-AD-vs-FirstReaction / Mark-calibration.")
        println(io)

        # --- matrix table ---
        println(io, "## Verdict Matrix (sampler x condition)")
        println(io)
        header = "| sampler | " * join(col_keys, " | ") * " |"
        sep = "|" * repeat("---|", length(col_keys) + 1)
        println(io, header)
        println(io, sep)
        for sk in sampler_order
            cells = String[]
            for ck in col_keys
                d = get(by_sc, (sk, ck), nothing)
                push!(cells, d === nothing ? "-" : _cell_triple(d))
            end
            println(io, "| ", sk, " | ", join(cells, " | "), " |")
        end
        println(io)

        # --- flagged cells ---
        println(io, "## Flagged Cells (BH-adjusted p <= $(alpha))")
        println(io)
        if isempty(flagged)
            println(io, "None. (See null-control discipline below for the milestone-1 raw flag.)")
        else
            println(io, "| sampler | condition | test | raw p | BH p | repro |")
            println(io, "|---|---|---|---|---|---|")
            for i in flagged
                r = runs[i]
                println(io, "| ", r.sampler_name, " | ", r.config_key, " | ",
                    TEST_SHORT[r.test_type], " | ", round(r.pvalue, digits=5), " | ",
                    round(adjusted[i], digits=5), " | `", repro_command(r; n=n), "` |")
            end
        end
        println(io)

        # --- raw sub-0.05 (pre-correction) for transparency ---
        raw_bugs = [(i, runs[i]) for i in eachindex(runs) if runs[i].verdict == :likely_bug]
        println(io, "## Raw sub-0.05 p-values (pre-correction, for transparency)")
        println(io)
        if isempty(raw_bugs)
            println(io, "None.")
        else
            println(io, "| sampler | condition | test | raw p | BH p | repro |")
            println(io, "|---|---|---|---|---|---|")
            for (i, r) in raw_bugs
                println(io, "| ", r.sampler_name, " | ", r.config_key, " | ",
                    TEST_SHORT[r.test_type], " | ", round(r.pvalue, digits=5), " | ",
                    isnan(adjusted[i]) ? "NA" : round(adjusted[i], digits=5),
                    " | `", repro_command(r; n=n), "` |")
            end
        end
        println(io)

        # --- grouping summaries ---
        println(io, "## Pattern Analysis (below-0.05 rate by parameter)")
        println(io)
        for param in (:test, :family, :memory, :graph, :sampler)
            g = group_by_parameter(runs, param)
            println(io, "### by `$param`")
            println(io)
            println(io, "| value | n | frac p<0.05 | min p |")
            println(io, "|---|---|---|---|")
            for k in sort(collect(keys(g)), by=string)
                s = g[k]
                println(io, "| ", k, " | ", s.n, " | ", round(s.frac_below_05, digits=3),
                    " | ", isnan(s.min_p) ? "NA" : round(s.min_p, digits=4), " |")
            end
            println(io)
        end

        # --- null-control investigations ---
        if !isempty(investigations)
            println(io, "## Null-Control Investigations")
            println(io)
            for inv in investigations
                println(io, "### ", inv.sampler, " @ ", inv.config)
                println(io)
                println(io, "Seed sweep over ", length(inv.seeds), " seeds. Fraction of seeds with p<0.05:")
                println(io)
                println(io, "| path | DM | AD | MK |")
                println(io, "|---|---|---|---|")
                println(io, "| SUT (", inv.sampler, ") | ", round(inv.sut_frac.doob, digits=2), " | ",
                    round(inv.sut_frac.two_sample, digits=2), " | ", round(inv.sut_frac.mark, digits=2), " |")
                println(io, "| null (FirstReaction vs itself) | ", round(inv.null_frac.doob, digits=2), " | ",
                    round(inv.null_frac.two_sample, digits=2), " | ", round(inv.null_frac.mark, digits=2), " |")
                println(io)
            end
        end

        println(io, "---")
        println(io, "Generated by `test/gauntlet/run_matrix.jl`.")
    end
    return path
end


"""
    print_summary_table(runs; alpha)

One-screen REPL summary: counts by band and the BH headline.
"""
function print_summary_table(runs::Vector{TestRun}; alpha::Float64=0.05)
    adjusted, flagged = apply_bh(runs; alpha=alpha)
    nbug = count(r -> r.verdict == :likely_bug, runs)
    nrerun = count(r -> r.verdict == :rerun, runs)
    nerr = count(r -> r.verdict == :error, runs)
    println("="^70)
    println("Matrix summary: $(length(runs)) results ($(length(runs) ÷ 3) cells x 3 tests)")
    println("  raw likely_bug (p<0.05): $nbug   rerun: $nrerun   errors: $nerr")
    println("  BH FDR (alpha=$alpha) flagged: $(length(flagged))")
    println("="^70)
    return nothing
end
