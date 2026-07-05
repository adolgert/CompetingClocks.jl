# matrix.jl - Phase 3 of the gauntlet: the full statistical test-matrix runner.
#
# Builds on the milestone-1 `run_gauntlet` (runner.jl), which already runs the
# three contract-conforming statistical tests for ONE sampler under ONE
# condition and returns three verdicts:
#   - :doob_meyer_stepcumulant          (Doob-Meyer step-cumulant uniformity)
#   - :two_sample_ad_vs_firstreaction   (contract-conforming two-sample vs FR)
#   - :mark_calibration                 (mark-calibration Brier permutation)
#
# The two-sample test uses `conditional_next_draws` (runner.jl) - the
# contract-conforming condensed-enabled-set path, NOT the out-of-contract
# replay path in experiments.jl - for ALL samplers.
#
# This file wraps `run_gauntlet` into a `TestConfiguration`-driven matrix over
# samplers x conditions, with a `TestRun` result record per (sampler, condition,
# test), deterministic per-cell seeds, Benjamini-Hochberg FDR correction across
# the whole matrix, and a null-control seed sweep for discipline.
#
# Assumes travel.jl, generate_data.jl, mark_calibration.jl, doob_meyer.jl,
# anderson_darling.jl, running_score.jl, runner.jl, and config.jl are in scope.

using Random


"""
    TestRun

One result of one statistical test in one matrix cell. `pvalue` is the raw
per-cell p-value; BH-adjusted p-values are computed across a whole run in
`benjamini_hochberg`. `verdict` is the plan's per-cell band. On an execution
error `pvalue` is `NaN`, `verdict` is `:error`, and `error_message` is set.
"""
struct TestRun
    sampler_name::String
    config_key::String
    config::TestConfiguration
    test_type::Symbol
    seed::UInt64
    pvalue::Float64
    verdict::Symbol
    n::Int
    elapsed::Float64
    error_message::Union{String,Nothing}
end


"Deterministic per-cell base seed from (base_seed, global cell index)."
cell_seed(base_seed::Integer, cell_index::Integer) = UInt64(base_seed) + UInt64(10 * cell_index)


"""
    run_single_cell(sampler_key, sampler_spec, config, n, seed) -> Vector{TestRun}

Run all three statistical tests for one (sampler, condition) cell by delegating
to `run_gauntlet` (contract-conforming path), relabeling the verdicts with the
matrix `sampler_key` (which distinguishes e.g. DirectMethod variants that
`run_gauntlet` would otherwise collapse to "DirectMethod"). Errors are captured
into `:error` TestRuns rather than aborting the matrix.
"""
function run_single_cell(sampler_key::AbstractString, sampler_spec, config::TestConfiguration,
        n::Int, seed::Integer)
    tcfg, builder = travel_config_and_builder(config)
    key = config_key(config)
    label = config_label(config)
    t0 = time()
    try
        verdicts = run_gauntlet(
            sampler_spec;
            n_replications=n,
            seed=seed,
            state_cnt=config.state_cnt,
            history_steps=config.history_steps,
            doob_steps=config.doob_steps,
            config=tcfg,
            model_builder=builder,
            condition_label=label,
            verbose=false,
        )
        elapsed = time() - t0
        return TestRun[
            TestRun(sampler_key, key, config, v.test, UInt64(seed),
                float(v.pvalue), v.verdict, v.n, elapsed, nothing)
            for v in verdicts
        ]
    catch err
        elapsed = time() - t0
        msg = sprint(showerror, err)
        return TestRun[
            TestRun(sampler_key, key, config, t, UInt64(seed),
                NaN, :error, 0, elapsed, msg)
            for t in (:doob_meyer_stepcumulant, :two_sample_ad_vs_firstreaction, :mark_calibration)
        ]
    end
end


"""
    run_matrix(; n, base_seed, samplers, design, verbose) -> Vector{TestRun}

Run the full sampler x condition matrix. For each sampler the design is filtered
to the conditions that sampler can run (`filter_for_sampler`); exponential-only
samplers (Direct*, RSSA, PSSACR) therefore run only the exponential family.

Cells are enumerated in a fixed order (registry order of samplers, then design
order of conditions), so the whole matrix is reproducible from `(base_seed, n)`.
Each cell gets a deterministic seed via `cell_seed`.
"""
function run_matrix(;
        n::Int=1000,
        base_seed::Integer=20260705,
        samplers=matrix_samplers(),
        design=generate_full_design(; doob_steps=n),
        verbose::Bool=true)
    runs = TestRun[]
    cell_index = 0
    ncells = sum(length(filter_for_sampler(design, cap)) for (_k, _s, cap) in samplers)
    verbose && println("Running matrix: $(ncells) cells x 3 tests, n=$(n), base_seed=$(base_seed)")
    for (skey, spec, cap) in samplers
        cfgs = filter_for_sampler(design, cap)
        for cfg in cfgs
            cs = cell_seed(base_seed, cell_index)
            cell_index += 1
            cell_runs = run_single_cell(skey, spec, cfg, n, cs)
            append!(runs, cell_runs)
            if verbose
                worst = minimum(r -> isnan(r.pvalue) ? Inf : r.pvalue, cell_runs)
                flag = any(r -> r.verdict == :likely_bug, cell_runs) ? " <-- flag" :
                       any(r -> r.verdict == :error, cell_runs) ? " <-- ERROR" : ""
                println(rpad("[$cell_index/$ncells] ", 12),
                    rpad(skey, 42), rpad(config_key(cfg), 26),
                    "min p=", rpad(round(worst, digits=4), 8), flag)
            end
        end
    end
    return runs
end


# --- Benjamini-Hochberg FDR correction --------------------------------------

"""
    benjamini_hochberg(pvals) -> Vector{Float64}

Benjamini-Hochberg step-up adjusted p-values (q-values) for a vector of raw
p-values. `NaN` entries (errored cells) are ignored in the ranking and returned
as `NaN`. adjusted[i] is comparable to a target FDR level: flag cells with
adjusted[i] <= alpha.
"""
function benjamini_hochberg(pvals::AbstractVector{<:Real})
    out = fill(NaN, length(pvals))
    idx = [i for i in eachindex(pvals) if !isnan(pvals[i])]
    m = length(idx)
    m == 0 && return out
    order = sort(idx, by=i -> pvals[i])
    # Step-up: adj_(k) = min_{j>=k} ( p_(j) * m / j ), then clamp to [0,1].
    running = 1.0
    for rank in m:-1:1
        i = order[rank]
        val = pvals[i] * m / rank
        running = min(running, val)
        out[i] = min(running, 1.0)
    end
    return out
end


"""
    apply_bh(runs; alpha) -> (adjusted::Vector{Float64}, flagged::Vector{Int})

Compute BH-adjusted p-values across ALL runs (the whole matrix is the family)
and return the indices of runs flagged at FDR level `alpha`.
"""
function apply_bh(runs::Vector{TestRun}; alpha::Float64=0.05)
    adjusted = benjamini_hochberg([r.pvalue for r in runs])
    flagged = [i for i in eachindex(runs) if !isnan(adjusted[i]) && adjusted[i] <= alpha]
    return adjusted, flagged
end


# --- Null-control discipline ------------------------------------------------

"""
    seed_sweep(sampler_key, config, seeds; n) -> Vector{NamedTuple}

Re-run one cell across several base seeds through the IDENTICAL path. Returns one
NamedTuple per seed with the three raw p-values. Used to decide whether a
BH-flag persists (a real finding) or was a single-seed fluctuation.
"""
function seed_sweep(sampler_key::AbstractString, config::TestConfiguration, seeds; n::Int=1000)
    spec = sampler_spec_for(sampler_key)
    out = NamedTuple[]
    for s in seeds
        runs = run_single_cell(sampler_key, spec, config, n, s)
        d = Dict(r.test_type => r.pvalue for r in runs)
        push!(out, (; seed=UInt64(s),
            doob=get(d, :doob_meyer_stepcumulant, NaN),
            two_sample=get(d, :two_sample_ad_vs_firstreaction, NaN),
            mark=get(d, :mark_calibration, NaN)))
    end
    return out
end


"""
    null_control(config, seeds; n) -> Vector{NamedTuple}

Run FirstReactionMethod as its OWN SUT (FirstReaction-vs-FirstReaction two-sample,
FirstReaction Doob-Meyer/mark) through the identical path under `config` across
`seeds`. This is the established null control: FR through the same machinery must
NOT flag more often than chance. Compare a flagged cell's seed sweep against this.
"""
function null_control(config::TestConfiguration, seeds; n::Int=1000)
    return seed_sweep("FirstReactionMethod", config, seeds; n=n)
end


"""
    investigate_flag(sampler_key, config; seeds, n) -> NamedTuple

Full null-control workup for a flagged cell: the flagged sampler's seed sweep and
the FirstReaction null-control sweep under the same condition, plus the fraction
of seeds at which each test dips below 0.05. A genuine finding shows the flagged
sampler dipping far more often than the FR null control.
"""
function investigate_flag(sampler_key::AbstractString, config::TestConfiguration;
        seeds=(20260705:20260705+9), n::Int=1000)
    sut = seed_sweep(sampler_key, config, seeds; n=n)
    null = null_control(config, seeds; n=n)
    frac_below(rows, field) = count(r -> !isnan(getproperty(r, field)) && getproperty(r, field) < 0.05, rows) / length(rows)
    return (;
        sampler=sampler_key, config=config_key(config),
        seeds=collect(seeds),
        sut, null,
        sut_frac=(; doob=frac_below(sut, :doob), two_sample=frac_below(sut, :two_sample), mark=frac_below(sut, :mark)),
        null_frac=(; doob=frac_below(null, :doob), two_sample=frac_below(null, :two_sample), mark=frac_below(null, :mark)),
    )
end
