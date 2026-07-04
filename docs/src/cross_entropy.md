# Cross-Entropy

Love it — here’s a clean **Cross-Entropy (CE) / Adaptive IS pattern** that plugs into your GSMP setup *without* needing analytic MLEs or per-clock sufficient stats. It treats your **bias multipliers** as decision variables, tests a small *population* of proposals each iteration, keeps the best (“elites”), and updates the proposal around those. It’s robust and easy to drop into your code.

---

# 0) What we’ll tune

Only the **first** element of each parameter pair is the *proposal*; the **second** stays the target. We’ll tune mild multipliers on a few knobs that don’t explode weights:

* `machine_off` (Exponential rate): ↓ rate ⇒ machine stays RUNNING longer.
* `fabricate_max` (burst size ceiling): ↑ slightly.
* `fabricate_setup` (machine setup rate): ↑ slightly.
* (Optionally) **avoid** tilting `assemble` (too many events); if you include it, keep changes tiny.

We’ll parameterize **log-multipliers** ( \phi=\log m ) (so (m=\exp \phi >0)) and do CE updates in log-space (numerically nicer, positivity guaranteed).

---

# 1) Minimal refactor: let `run_epochs` accept a param dictionary

Your `run_epochs` currently bakes the params inside. Make it accept a `params` Dict so the CE loop can swap them per candidate:

```julia
# --- small refactor ---
function run_epochs(epoch_cnt::Int, importance::Bool, rng::AbstractRNG, params::Dict)
    if !importance
        # erase weighted params (proposal = target) for plain simulation
        params = Dict((k => (v[2], v[2])) for (k, v) in params)
    end
    widget = zeros(Int, epoch_cnt)
    importance = zeros(Float64, epoch_cnt)
    for epoch_idx in eachindex(widget)
        model = Workstation(params)
        builder = SamplerBuilder(
            Tuple{Symbol,Int}, Float64;
            sampler_spec=:firsttofire,
            trajectory_likelihood=true,
            likelihood_cnt=2,
        )
        sampler = SamplingContext(builder, rng)
        (cnt, w) = one_epoch(model, sampler)   # w = exp(basal - weighted)
        widget[epoch_idx] = cnt
        importance[epoch_idx] = w
    end
    return widget, importance
end
```

And keep a **baseline** target set somewhere central:

```julia
const TARGET_PARAMS = Dict(
    :machine_off         => (missing, 0.6),   # per min  (proposal gets filled)
    :fabricate_max       => (missing, 10.0),  # parts/min
    :fabricate_setup     => (missing, 1.0),   # per min
    :degrade_k           => (missing, 4.0),
    :degrade_theta       => (missing, 4 * 4 / 30),
    :assemble            => (missing, 1.0)    # widgets/min/part
)
```

> Note: only the **second** entry is used as the target. We’ll construct the first entry (proposal) from multipliers.

---

# 2) Utility: build a proposal from log-multipliers

```julia
"""
make_params(φ; include_assemble=false)

φ is a NamedTuple of log-multipliers, e.g.
    (machine_off = log(0.4/0.6), fabricate_max = log(1.2), fabricate_setup = log(1.5), assemble = log(1.0))

Returns a Dict matching the (proposal, target) shape your model expects.
"""
function make_params(φ; include_assemble::Bool=false)
    tgt = TARGET_PARAMS
    # Start proposal equal to target
    prop = Dict(
        :machine_off         => tgt[:machine_off][2],
        :fabricate_max       => tgt[:fabricate_max][2],
        :fabricate_setup     => tgt[:fabricate_setup][2],
        :degrade_k           => tgt[:degrade_k][2],
        :degrade_theta       => tgt[:degrade_theta][2],
        :assemble            => tgt[:assemble][2],
    )
    # Apply safe tilts (multiplicative in *rate* space)
    prop[:machine_off]        *= exp(φ.machine_off)
    prop[:fabricate_max]      *= exp(φ.fabricate_max)
    prop[:fabricate_setup]    *= exp(φ.fabricate_setup)
    if include_assemble
        prop[:assemble]      *= exp(φ.assemble)
    end
    return Dict(
        :machine_off         => (prop[:machine_off],        tgt[:machine_off][2]),
        :fabricate_max       => (prop[:fabricate_max],      tgt[:fabricate_max][2]),
        :fabricate_setup     => (prop[:fabricate_setup],    tgt[:fabricate_setup][2]),
        :degrade_k           => (prop[:degrade_k],          tgt[:degrade_k][2]),
        :degrade_theta       => (prop[:degrade_theta],      tgt[:degrade_theta][2]),
        :assemble            => (prop[:assemble],           tgt[:assemble][2]),
    )
end
```

---

# 3) Scoring a candidate proposal

We’ll run (N) paths, use **log-space** + **self-normalized IS**, and compute a score that trades off **success frequency** and **weight quality**. A simple, effective scalar score is:

[
\text{score} ;=; \text{ESS}*\text{norm} \times \text{fraction_over}
]
where (\text{ESS}*\text{norm}=\tfrac{(\sum w)^2}{N\sum w^2}\in[0,1]).

```julia
using Random, Statistics

struct CandidateResult
    φ::NamedTuple
    frac_over::Float64
    p_over::Float64   # SNIS estimate
    ess_norm::Float64 # in [0,1]
    score::Float64
end

function evaluate_candidate(φ; N=5_000, rng=Xoshiro(0), threshold=1000, include_assemble=false)
    params = make_params(φ; include_assemble)
    obs, Δ = run_epochs(N, true, rng, params)
    # Stable SNIS weights
    mΔ = maximum(Δ); w = exp.(Δ .- mΔ)
    frac_over = count(>(threshold), obs) / N
    p_over = sum((obs .>= threshold) .* w) / sum(w)
    ess_norm = (sum(w)^2) / (N * sum(w.^2))
    score = ess_norm * frac_over
    return CandidateResult(φ, frac_over, p_over, ess_norm, score)
end
```

---

# 4) CE loop over **log-multipliers**

We keep a Gaussian on φ (log-multipliers), sample a small population, pick the top ρ-quantile, and update the mean (and optionally the std) with **smoothing** to avoid oscillations.

```julia
"""
cross_entropy_opt

- φμ, φσ: NamedTuples of means/stds in LOG space for the multipliers to tune
- K:      population size per iteration
- ρ:      elite fraction (e.g., 0.2)
- α:      smoothing for mean updates (0<α≤1)
- target_frac: optional target range for fraction_over (keeps proposals "gentle")

Returns the best proposal found and a log of iterations.
"""
function cross_entropy_opt(; 
    φμ = (machine_off=log(0.2/0.6), fabricate_max=log(1.0), fabricate_setup=log(1.0), assemble=log(1.0)),
    φσ = (machine_off=0.3,           fabricate_max=0.15,    fabricate_setup=0.15,     assemble=0.0),
    include_assemble=false,
    iters=8, K=12, ρ=0.25, α=0.6, N=5_000, rng=Xoshiro(123),
    target_frac=(0.05, 0.15), threshold=1000
)
    results_log = Vector{CandidateResult}()
    best::Union{Nothing,CandidateResult} = nothing

    for t in 1:iters
        cand = Vector{CandidateResult}(undef, K)
        for k in 1:K
            φ = (
                machine_off     = randn(rng) * φσ.machine_off     + φμ.machine_off,
                fabricate_max   = randn(rng) * φσ.fabricate_max   + φμ.fabricate_max,
                fabricate_setup = randn(rng) * φσ.fabricate_setup + φμ.fabricate_setup,
                assemble        = include_assemble ? (randn(rng) * φσ.assemble + φμ.assemble) : 0.0,
            )
            cand[k] = evaluate_candidate(φ; N=N, rng=rng, threshold=threshold, include_assemble=include_assemble)
        end
        # Optionally filter to keep fraction_over in a sane window
        cand_filtered = filter(c -> target_frac[1] <= c.frac_over <= target_frac[2], cand)
        pool = isempty(cand_filtered) ? cand : cand_filtered

        # Pick elites by score
        sort!(pool, by = c -> c.score, rev=true)
        elites = pool[1:clamp(ceil(Int, ρ*length(pool)), 1, length(pool))]

        # Update φμ via elite mean (optionally weight by score)
        new_μ = (
            machine_off     = mean(e.φ.machine_off     for e in elites),
            fabricate_max   = mean(e.φ.fabricate_max   for e in elites),
            fabricate_setup = mean(e.φ.fabricate_setup for e in elites),
            assemble        = include_assemble ? mean(e.φ.assemble for e in elites) : 0.0,
        )
        # Smooth update
        φμ = (
            machine_off     = (1-α)*φμ.machine_off     + α*new_μ.machine_off,
            fabricate_max   = (1-α)*φμ.fabricate_max   + α*new_μ.fabricate_max,
            fabricate_setup = (1-α)*φμ.fabricate_setup + α*new_μ.fabricate_setup,
            assemble        = include_assemble ? ((1-α)*φμ.assemble + α*new_μ.assemble) : 0.0,
        )

        # (Optional) shrink stds a bit to concentrate
        φσ = (
            machine_off     = 0.9 * φσ.machine_off,
            fabricate_max   = 0.9 * φσ.fabricate_max,
            fabricate_setup = 0.9 * φσ.fabricate_setup,
            assemble        = include_assemble ? (0.9 * φσ.assemble) : 0.0,
        )

        # Track the single best seen
        for e in cand
            if best === nothing || e.score > best.score
                best = e
            end
        end

        # Log a line for this iteration (best in this gen)
        println("-- iter $t  best score=$(pool[1].score): frac_over=$(round(pool[1].frac_over, digits=4))  p_over=$(round(pool[1].p_over, sigdigits=3))  ESS%=$(round(100*pool[1].ess_norm, digits=1))")
        append!(results_log, cand)
    end

    @assert best !== nothing
    return best, results_log, φμ, φσ
end
```

**Usage:**

```julia
best, logres, φμ_final, φσ_final = cross_entropy_opt(
    φμ=(machine_off=log(0.2/0.6), fabricate_max=log(1.0), fabricate_setup=log(1.0), assemble=log(1.0)),
    φσ=(machine_off=0.25,          fabricate_max=0.15,      fabricate_setup=0.15,     assemble=0.0),
    include_assemble=false,   # safer
    iters=8, K=12, ρ=0.25, α=0.6, N=5_000, rng=Xoshiro(42),
    target_frac=(0.05, 0.15), threshold=1000
)

println("\nBest proposal found:")
println("  frac_over = $(best.frac_over)")
println("  p_over(SNIS) = $(best.p_over)")
println("  ESS% = $(round(100*best.ess_norm, digits=1))")
println("  log-multipliers φ = ", best.φ)
println("  multipliers m = ", (machine_off=exp(best.φ.machine_off), fabricate_max=exp(best.φ.fabricate_max), fabricate_setup=exp(best.φ.fabricate_setup)))
```

---

## Why this pattern works well for your GSMP

* It **never** touches your inner likelihood recorder — we keep using your correct `PathLikelihoods` to compute Δ for both target and proposal.
* It **avoids** the hard analytic refits (which would need per-clock sufficient stats) by using a **black-box** objective (ESS × success frequency).
* It **prefers gentle proposals** by (a) filtering to a target success window (e.g., 5–15%), and (b) shrinking stds each iteration.
* It’s **numerically stable** (log-space weights, SNIS).
* It’s **modular**: add/remove bias knobs by extending `φμ`, `φσ`, and `make_params`.

---

## A couple of practical tips

* Start with **machine_off** and **fabricate_setup**; leave **assemble** off (variance bomb).
* Keep **N modest** (3–10k) during CE search; when you lock a proposal, do one **large** run to finalize the estimate and confidence intervals.
* Watch **ESS%**; you generally want ≥30–40% while maintaining a healthy `frac_over` (5–15% is a nice target).

If you want, we can tailor the knobs further (e.g., only bias **early** machine warm-up by using a shifted distribution for `SetupRate` in the first X minutes). That tends to improve efficiency even more for “early burst” rare events.

