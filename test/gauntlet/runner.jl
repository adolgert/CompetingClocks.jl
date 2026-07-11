# runner.jl - Milestone 1 of the statistical-correctness gauntlet.
#
# This runs the gauntlet END TO END for ONE sampler under ONE condition (the
# simplest one: all-Exponential rates, forget/no-memory, no enabling delay) and
# produces statistical verdicts. It also supports MultiSampler compositions
# (sub-samplers keyed by a chooser predicate) as the sampler under test, and a
# pluggable `model_builder` so a non-exponential (Weibull) TravelModel variant
# can be tested without modifying travel.jl.
#
# It reuses the package's own machinery and the existing gauntlet building
# blocks rather than reimplementing statistics:
#   - Doob-Meyer / step-cumulant uniformity test: driven by the `ScoreEvents`
#     observer (running_score.jl), which calls `CompetingClocks.stepcumulant`
#     on a `CompetingClocks.TrackWatcher` recording context. The collected
#     step-cumulants should be Uniform(0,1); uniformity is checked with
#     `HypothesisTests.ApproximateOneSampleKSTest`.
#   - Two-sample test vs the trusted FirstReactionMethod reference: condenses a
#     fixed command history into its final enabled set (via the existing
#     `final_enabled_distributions` helper), draws one next-event from each of
#     many fresh replicas of both samplers, and compares the firing-time
#     distributions with the two-sample Anderson-Darling test
#     `HypothesisTests.KSampleADTest`. See `conditional_next_draws` for why the
#     query is structured this way (the `next` contract).
#
# This file assumes the sibling gauntlet files have already been included in the
# enclosing scope (travel.jl, generate_data.jl, mark_calibration.jl,
# running_score.jl, doob_meyer.jl, anderson_darling.jl). See
# test_gauntlet_smoke.jl for a self-contained entry point.

using CompetingClocks
using CompetingClocks: NextReactionMethod, FirstReactionMethod, FirstToFireMethod
using HypothesisTests
using Distributions
using Random


"""
    simple_condition()

The single, simplest gauntlet condition for milestone 1: exponentially
distributed rates, no memory (forget), no enabling delay/shift, one rate per
destination, on a cycle graph. This matches the plan's "all-Exponential,
no memory/forgetting" baseline.
"""
function simple_condition()
    return TravelConfig(
        TravelMemory.forget,
        TravelGraph.cycle,
        TravelRateDist.exponential,
        TravelRateCount.destination,
        TravelRateDelay.none,
    )
end


"""
    pvalue_verdict(p) -> Symbol

Map a p-value to the plan's interpretation bands:
  - `p < 0.05`  -> `:likely_bug`      (statistically significant discrepancy)
  - `p > 0.10`  -> `:likely_correct`
  - otherwise   -> `:rerun`           (marginal; rerun with new randomness)
"""
function pvalue_verdict(p::Real)
    if p < 0.05
        return :likely_bug
    elseif p > 0.10
        return :likely_correct
    else
        return :rerun
    end
end


"""
    multisampler_condition()

Condition for MultiSampler-under-test runs: same as `simple_condition()` but on
a COMPLETE graph. With a complete graph on `state_cnt >= 4` sites the enabled
set at every step is all destinations except the current site, so a parity
split of the clocks always places clocks in BOTH sub-sampler groups — including
in the condensed enabled set used by the two-sample test. (On the 5-cycle, the
2-clock enabled sets can be single-parity, which would leave one group idle.)
"""
function multisampler_condition()
    return TravelConfig(
        TravelMemory.forget,
        TravelGraph.complete,
        TravelRateDist.exponential,
        TravelRateCount.destination,
        TravelRateDelay.none,
    )
end


"""
    default_model_builder(state_cnt, config, rng) -> Travel

The standard TravelModel construction (exponential per-destination rates, as
`Travel(state_cnt, config, rng)` builds today).
"""
default_model_builder(state_cnt::Int, config::TravelConfig, rng::AbstractRNG) =
    Travel(state_cnt, config, rng)


"""
    weibull_model_builder(state_cnt, config, rng) -> Travel

A non-exponential TravelModel variant: one Weibull distribution per destination
(shapes spread across [0.7, 1.6], so both decreasing- and increasing-hazard
clocks appear), no enabling delay. Built directly from the `Travel` struct
constructor because the existing `Travel(state_cnt, config, rng)` constructor
ignores `config.dist` and always builds exponential rates (see travel.jl);
wiring `TravelRateDist.general` through that constructor is matrix-milestone
work. Non-exponential rates are where the two-sample test has real power:
inside the contract, CombinedNextReaction's inversion for non-exponentials is
not the same uniform-to-time transform as FirstReaction's, so the comparison is
no longer trivially pathwise-identical.
"""
function weibull_model_builder(state_cnt::Int, config::TravelConfig, rng::AbstractRNG)
    g = travel_make_graph(config.graph, state_cnt)
    shapes = collect(range(0.7, stop=1.6, length=state_cnt))
    scales = exp.(range(-1.0, stop=1.0, length=state_cnt))
    rates = Dict{Tuple{Int,Int},Tuple{UnivariateDistribution,Float64}}()
    for i in 1:state_cnt, j in 1:state_cnt
        i == j && continue
        rates[(i, j)] = (Weibull(shapes[j], scales[j]), 0.0)
    end
    return Travel(g, rates, zeros(Float64, state_cnt), config.memory)
end


# --- MultiSampler-under-test specs -----------------------------------------
#
# The runner builds samplers by calling `sampler_spec(Int, Float64)`. These
# specs extend that interface to CompetingClocks.MultiSampler compositions:
# a chooser (SamplerChoice) routes each clock to a sub-sampler group by a key
# predicate, and each group holds its own low-level sampler.

"Route travel clocks (destination site numbers) to `:even`/`:odd` groups."
struct GauntletParitySplit <: CompetingClocks.SamplerChoice{Symbol,Int} end
function CompetingClocks.choose_sampler(
    ::GauntletParitySplit, clock::Int, ::UnivariateDistribution
)::Symbol
    return iseven(clock) ? :even : :odd
end

"Route every clock to the single `:all` group (control composition)."
struct GauntletOneGroup <: CompetingClocks.SamplerChoice{Symbol,Int} end
function CompetingClocks.choose_sampler(
    ::GauntletOneGroup, clock::Int, ::UnivariateDistribution
)::Symbol
    return :all
end


"""
    MultiSamplerGauntletSpec(name, chooser, groups)

A sampler-spec-like callable for the gauntlet runner: `spec(K, T)` builds a
`CompetingClocks.MultiSampler{Symbol,K,T}` with the given chooser and one
sub-sampler per `group_key => sub_spec` pair (each `sub_spec` is itself a
sampler spec such as `NextReactionMethod()`).
"""
struct MultiSamplerGauntletSpec{C}
    name::String
    chooser::C
    groups::Vector
end

function (spec::MultiSamplerGauntletSpec)(K, T)
    ms = CompetingClocks.MultiSampler{Symbol,K,T}(spec.chooser)
    for (group_key, sub_spec) in spec.groups
        ms[group_key] = sub_spec(K, T)
    end
    return ms
end

"Condition (a): two groups split by destination parity, both CombinedNextReaction."
multisampler_split_spec() = MultiSamplerGauntletSpec(
    "MultiSampler(even=>CNR,odd=>CNR)",
    GauntletParitySplit(),
    [:even => NextReactionMethod(), :odd => NextReactionMethod()],
)

"Condition (b): mixed sub-samplers, CombinedNextReaction and FirstToFire."
multisampler_mixed_spec() = MultiSamplerGauntletSpec(
    "MultiSampler(even=>CNR,odd=>FirstToFire)",
    GauntletParitySplit(),
    [:even => NextReactionMethod(), :odd => FirstToFireMethod()],
)

"Condition (c): control, ONE group containing all clocks (CombinedNextReaction)."
multisampler_control_spec() = MultiSamplerGauntletSpec(
    "MultiSampler(all=>CNR)",
    GauntletOneGroup(),
    [:all => NextReactionMethod()],
)


"Human-readable sampler name for the verdict table."
gauntlet_sampler_name(spec) = string(nameof(typeof(spec)))
gauntlet_sampler_name(spec::MultiSamplerGauntletSpec) = spec.name


"""
    doob_meyer_stepcumulant(sampler_spec, config, state_cnt, n_steps, rng)

Run a single trajectory of `n_steps` steps of the sampler under test and collect
the package's step-cumulants via `ScoreEvents` (which calls
`CompetingClocks.stepcumulant` on a `CompetingClocks.TrackWatcher`). Each
step-cumulant is Uniform(0,1) under a correct sampler. All per-state cumulants
are pooled and tested for uniformity with a one-sample KS test. The
`ScoreEvents` observer also produces a mark-calibration score in the same pass,
which we surface as a bonus.

Returns a NamedTuple with `pvalue` (pooled Doob-Meyer), `n`, per-state detail,
and the mark-calibration result.
"""
function doob_meyer_stepcumulant(
    sampler_spec, config, state_cnt::Int, n_steps::Int, rng::Xoshiro;
    model_builder = default_model_builder,
)
    sampler = sampler_spec(Int, Float64)
    # Seed the sampler's own streams from the passed rng so this trajectory is
    # a genuine random draw rather than the default-seed one.
    CompetingClocks.rekey_streams!(sampler, rand(rng, UInt64))
    model = model_builder(state_cnt, config, rng)
    tracker = CompetingClocks.TrackWatcher{Int64,Float64}()
    observer = ScoreEvents(tracker, rng)
    # SamplerRecord mirrors every enable!/disable!/fire! into the TrackWatcher so
    # the observer sees the live set of enabled clocks (the recording context).
    recorder = SamplerRecord(tracker, rng)
    travel_run(n_steps, sampler, model, observer, recorder, rng)

    # Pool all step-cumulants across originating states. Each is Uniform(0,1)
    # marginally under a correct sampler, so the pool should be Uniform(0,1).
    pooled = Float64[]
    for v in values(observer.waiting_cumulant)
        append!(pooled, v)
    end
    test = ApproximateOneSampleKSTest(pooled, Uniform(0, 1))
    p = pvalue(test)

    per_state = waiting_metric(observer)   # per-state Doob-Meyer p-values
    mark = mark_calibration(observer)      # bonus: mark-calibration score
    return (; pvalue=p, n=length(pooled), supremum=test.δ, per_state, mark)
end


"""
    conditional_next_draws(sampler_spec, enabled, t_final, n_replications, rng)

Draw `n_replications` samples of the next event `(clock, time)` given a history
whose final event was at `t_final` and whose surviving clocks are `enabled`
(a `Dict{Int,DistributionState}` mapping clock to distribution and original
enabling time `te <= t_final`).

CONTRACT NOTE (see the `next` docstring in `src/sample/interface.jl`): the
`when` argument of `next(sampler, when)` must never decrease and must
never advance past a pending firing time without that event being fired; the
time of the most recently fired event is always a safe choice. Outside that
invariant samplers legitimately disagree (FirstReaction re-conditions on
survival at query time; CombinedNextReaction returns putative times cached at
`enable!` time). Replaying a recorded history into a fresh replica sampler and
then querying at the history's final time violates the invariant: the replica
draws its own putative times at each `enable!`, and the replayed history
advances the clock past them without firing them. So instead we condense the
history into its final enabled set and, for each replica, call
`enable!(sampler, clock, dist, te, t_final, rng)` with the ORIGINAL enabling
time `te` (possibly in the past) at the CURRENT time `t_final`, then query
`next(sampler, t_final)`. Now `t_final` is exactly the time of the most
recent state change, no pending event predates it, and conditioning on survival
through the history is handled uniformly inside `enable!` (both samplers
truncate to `t_final - te` when `te < t_final`). This measures what the plan's
two-sample test claims to measure: the distribution of the next event given the
history.
"""
function conditional_next_draws(
    sampler_spec, enabled::AbstractDict, t_final::Float64,
    n_replications::Int, rng::AbstractRNG,
)
    draws = Vector{ClockDraw}(undef, n_replications)
    for rep in 1:n_replications
        sampler = sampler_spec(Int, Float64)
        # The sampler owns its randomness; give each replication its own seed
        # (drawn from the passed rng) so the reps are independent draws of the
        # conditional next-event rather than identical repeats.
        CompetingClocks.rekey_streams!(sampler, rand(rng, UInt64))
        for (clock, ds) in enabled
            enable!(sampler, clock, ds.d, ds.enabling_time, t_final)
        end
        when, which = next(sampler, t_final)
        @assert which !== nothing
        @assert when >= t_final
        draws[rep] = (which, when)
    end
    return draws
end


"""
    gauntlet_history_commands(model_builder, history_steps, state_cnt, config, rng)

Generate a fixed command history: a `history_steps`-step FirstReaction
trajectory of the model built by `model_builder`. With the default builder this
matches `TravelModel.travel_commands`; the builder injection lets the same
history generation serve non-exponential model variants.
"""
function gauntlet_history_commands(model_builder, history_steps::Int, state_cnt::Int, config, rng)
    model = model_builder(state_cnt, config, rng)
    sampler = CompetingClocks.FirstReaction{Int,Float64}()
    rec = VectorRecord()
    travel_run(history_steps, sampler, model, (x, y) -> nothing, rec, rng)
    return rec.commands
end


"""
    two_sample_ad_vs_reference(sampler_spec, config, state_cnt, history_steps,
                               n_replications, rng)

Build a fixed command history (a `history_steps`-step FirstReaction trajectory),
condense it to its final enabled set, then draw `n_replications` next-events
from both the trusted `FirstReactionMethod` reference and the sampler under
test with `conditional_next_draws` (which documents why the query is structured
to conform to the `next` contract). The next-event firing-time distributions
are compared with the two-sample Anderson-Darling test
(`HypothesisTests.KSampleADTest`), pooled over all clocks and also per clock.

Returns a NamedTuple with the aggregate `pvalue`, `n`, and per-clock detail.
"""
function two_sample_ad_vs_reference(
    sampler_spec, config, state_cnt::Int, history_steps::Int,
    n_replications::Int, rng::AbstractRNG;
    model_builder = default_model_builder,
)
    commands = gauntlet_history_commands(model_builder, history_steps, state_cnt, config, rng)
    enabled = final_enabled_distributions(commands)
    t_final = commands[end][end]::Float64

    draws_ref = conditional_next_draws(FirstReactionMethod(), enabled, t_final, n_replications, rng)
    draws_sut = conditional_next_draws(sampler_spec, enabled, t_final, n_replications, rng)

    # Aggregate over all clocks: compare the marginal next-firing-time samples.
    times_ref = Float64[x[2] for x in draws_ref]
    times_sut = Float64[x[2] for x in draws_sut]
    agg = KSampleADTest(times_ref, times_sut)
    p_agg = pvalue(agg)

    # Per-clock comparison of conditional firing times.
    clocks = sort(unique(x[1] for x in vcat(draws_ref, draws_sut)))
    per_clock = NamedTuple[]
    for c in clocks
        ta = Float64[x[2] for x in draws_ref if x[1] == c]
        tb = Float64[x[2] for x in draws_sut if x[1] == c]
        (length(ta) < 2 || length(tb) < 2) && continue
        t = KSampleADTest(ta, tb)
        push!(per_clock, (; clock=c, pvalue=pvalue(t), n_ref=length(ta), n_sut=length(tb)))
    end

    return (; pvalue=p_agg, n=n_replications, per_clock)
end


"""
    run_gauntlet(sampler_spec=NextReactionMethod(); n_replications, seed, ...)

Run the milestone-1 gauntlet for one `sampler_spec` under the single simplest
condition and return a `Vector` of verdict NamedTuples, one per statistical test.
Each verdict has fields `(sampler, test, pvalue, verdict, n, seed)`.

Keyword arguments:
  - `n_replications`: primary sample-size knob. Drives both the Doob-Meyer
    trajectory length (`doob_steps`) and the number of two-sample replicas.
  - `seed`: base RNG seed; all randomness derives deterministically from it.
  - `state_cnt`: number of sites in the TravelModel (default 5).
  - `history_steps`: length of the fixed history for the two-sample test.
  - `doob_steps`: Doob-Meyer trajectory length (defaults to `n_replications`).
  - `config`: the TravelModel condition (defaults to `simple_condition()`).
  - `model_builder`: how to construct the Travel model from `(state_cnt,
    config, rng)`; pass `weibull_model_builder` for the non-exponential variant.
  - `condition_label`: condition description for the verdict table header.
  - `verbose`: print the verdict table to stdout.
"""
function run_gauntlet(
    sampler_spec = NextReactionMethod();
    n_replications::Int = 1000,
    seed::Integer = 0x5eed,
    state_cnt::Int = 5,
    history_steps::Int = 5,
    doob_steps::Int = n_replications,
    config::TravelConfig = simple_condition(),
    model_builder = default_model_builder,
    condition_label::AbstractString = string(
        "memory=", config.memory, " graph=", config.graph,
        " delay=", config.delay, " rates=",
        model_builder === weibull_model_builder ? "Weibull-per-destination" : "Exponential-per-destination",
    ),
    verbose::Bool = true,
)
    sampler_name = gauntlet_sampler_name(sampler_spec)

    dm = doob_meyer_stepcumulant(
        sampler_spec, config, state_cnt, doob_steps, Xoshiro(seed);
        model_builder,
    )
    ts = two_sample_ad_vs_reference(
        sampler_spec, config, state_cnt, history_steps, n_replications,
        Xoshiro(seed + 1);
        model_builder,
    )

    verdicts = [
        (; sampler=sampler_name, test=:doob_meyer_stepcumulant,
            pvalue=dm.pvalue, verdict=pvalue_verdict(dm.pvalue), n=dm.n, seed=UInt64(seed)),
        (; sampler=sampler_name, test=:two_sample_ad_vs_firstreaction,
            pvalue=ts.pvalue, verdict=pvalue_verdict(ts.pvalue), n=ts.n, seed=UInt64(seed)),
        (; sampler=sampler_name, test=:mark_calibration,
            pvalue=dm.mark.pvalue, verdict=pvalue_verdict(dm.mark.pvalue), n=dm.mark.count, seed=UInt64(seed)),
    ]

    if verbose
        print_verdict_table(verdicts; sampler_name, seed, n_replications, dm, ts, condition_label)
    end
    return verdicts
end


"""
    print_verdict_table(verdicts; ...)

Print a small, human-readable verdict table to stdout with the plan's
interpretation bands.
"""
function print_verdict_table(
    verdicts;
    sampler_name="", seed=0, n_replications=0, dm=nothing, ts=nothing,
    condition_label="Exponential rates, memory=forget, delay=none, graph=cycle, per-destination",
)
    println("=" ^ 78)
    println("Gauntlet verdicts - sampler=$(sampler_name)")
    println("condition: $(condition_label)")
    println("n_replications=$(n_replications)  seed=$(seed)")
    println("-" ^ 78)
    println(rpad("test", 34), lpad("p-value", 10), "  ", lpad("n", 6), "  ", "verdict")
    println("-" ^ 78)
    for v in verdicts
        println(
            rpad(string(v.test), 34),
            lpad(string(round(v.pvalue, digits=4)), 10), "  ",
            lpad(string(v.n), 6), "  ",
            string(v.verdict),
        )
    end
    println("-" ^ 78)
    println("bands: p<0.05 => likely_bug   0.05<=p<=0.10 => rerun   p>0.10 => likely_correct")
    if dm !== nothing && !isempty(dm.per_state)
        worst = minimum(s.pvalue for s in dm.per_state)
        println("  doob per-state: $(length(dm.per_state)) states, min p-value=$(round(worst, digits=4))")
    end
    if ts !== nothing && !isempty(ts.per_clock)
        worst = minimum(s.pvalue for s in ts.per_clock)
        println("  two-sample per-clock: $(length(ts.per_clock)) clocks, min p-value=$(round(worst, digits=4))")
    end
    println("=" ^ 78)
    return nothing
end
