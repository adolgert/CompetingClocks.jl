# Tests for the estimator-facing sampler verbs (Milestone 3): the capability
# traits, `enabled_ages`, `retained_draw`, and `force_fire!`. Plain @testset
# (not @safetestset) so the file-scope helpers and quadrature oracles below are
# in scope for every testset, matching the mixed style of test_recorder.jl.
#   julia --project=test test/runtests.jl "estimator_verbs"
using Test
using CompetingClocks
using CompetingClocks: FirstReaction, FirstToFire, CombinedNextReaction,
    DirectCall, MultipleDirect, RSSA, PSSACR, Petri,
    supports_enabled_ages, supports_force, supports_retained_draw,
    enabling_times, enabled_ages, retained_draw, force_fire!,
    enable!, disable!, fire!, next, isenabled
using Distributions: Weibull, Exponential, Uniform, LogNormal, Normal, Gamma,
    truncated, logpdf, logccdf, invlogccdf, ccdf, pdf
using QuadGK: quadgk
using Random: Xoshiro
using Statistics: mean, std


# --- quadrature oracles for the forced-firing conditional law ---------------
# The exact law the survivors of a forced firing at `tstar` must follow: each
# survivor's lifetime, measured from its enabling time `te`, conditioned on
# survival past `tstar`. `clocks` is an indexable collection of (dist, te)
# survivors. Ported from src/ForcedFire/oracles.jl so the law test compares both
# backends against ground truth, not only against each other.

ev_survivor_logccdf(dist, te, tstar, t) =
    logccdf(dist, t - te) - logccdf(dist, tstar - te)

# P(survivor i fires first after the force) =
#   ∫_{t*}^∞ [pdf_i(t−te_i)/S_i(t*−te_i)] ∏_{j≠i} S̃_j(t) dt.
function ev_win_probability(clocks, i, tstar; rtol=1e-9)
    di, tei = clocks[i].dist, clocks[i].te
    quadgk(Float64(tstar), Inf; rtol=rtol) do t
        lw = logpdf(di, t - tei) - logccdf(di, tstar - tei)
        for (j, c) in enumerate(clocks)
            j == i && continue
            lw += ev_survivor_logccdf(c.dist, c.te, tstar, t)
        end
        exp(lw)
    end[1]
end

# E[next firing time] = t* + ∫_{t*}^∞ ∏_k S̃_k(t) dt (tail-integral identity).
function ev_mean_time(clocks, tstar; rtol=1e-9)
    Float64(tstar) + quadgk(Float64(tstar), Inf; rtol=rtol) do t
        exp(sum(ev_survivor_logccdf(c.dist, c.te, tstar, t) for c in clocks))
    end[1]
end

# E[T · 1{i wins next}] — a joint (winner, time) moment.
function ev_joint_moment(clocks, i, tstar; rtol=1e-9)
    di, tei = clocks[i].dist, clocks[i].te
    quadgk(Float64(tstar), Inf; rtol=rtol) do t
        lw = logpdf(di, t - tei) - logccdf(di, tstar - tei)
        for (j, c) in enumerate(clocks)
            j == i && continue
            lw += ev_survivor_logccdf(c.dist, c.te, tstar, t)
        end
        t * exp(lw)
    end[1]
end

# The three backends that carry a per-clock enabling-time table and therefore
# support the estimator verbs. FirstReaction is conditional-redraw; the other
# two are scheduling backends.
const EV_AGED_BACKENDS = (
    ("FirstReaction", FirstReaction),
    ("FirstToFire", FirstToFire),
    ("CombinedNextReaction", CombinedNextReaction),
)


@testset "estimator_verbs: enabled_ages returns every enabled clock with its age when - te, sorted by key, including a left-shifted enable whose age exceeds the elapsed simulation time, on all three aged backends" begin
    # Enable at simulation time when = 1.0. Clock :z is left-shifted (te = -0.4),
    # so its age 1.4 exceeds the 1.0 of elapsed simulation time — the case an
    # estimator must get right because a left-shifted clock is aged before the
    # sampler ever saw it.
    when = 1.0
    dists = Dict(:a => Weibull(1.7, 1.0), :m => Exponential(1.2),
                 :z => Weibull(2.0, 1.5))
    tes = Dict(:a => 0.0, :m => 0.5, :z => -0.4)
    expected = sort([(k, when - tes[k]) for k in keys(tes)]; by=first)

    for (label, B) in EV_AGED_BACKENDS
        rng = Xoshiro(4040)
        s = B{Symbol,Float64}()
        for k in (:z, :a, :m)   # deliberately not in sorted order
            enable!(s, k, dists[k], tes[k], when)
        end
        ages = enabled_ages(s, when)
        # Sorted by key regardless of enable order (deterministic order is what
        # forced-selection pmfs index into).
        @test [k for (k, _) in ages] == [k for (k, _) in expected]
        for ((k, age), (ek, eage)) in zip(ages, expected)
            @test k == ek
            @test isapprox(age, eage; atol=1e-12)
        end
        @test any(age > when for (_, age) in ages)   # the left-shifted clock
    end
end


@testset "estimator_verbs: enabled_ages excludes fired and disabled clocks, and on CombinedNextReaction the retained-survival entries kept for draw reuse (heap_handle == 0) do not reappear" begin
    for (label, B) in EV_AGED_BACKENDS
        rng = Xoshiro(5050)
        s = B{Symbol,Float64}()
        for k in (:a, :b, :c, :d)
            enable!(s, k, Exponential(1.0), 0.0, 0.0)
        end
        disable!(s, :b, 0.2)      # cancelled
        fire!(s, :a, 0.3)         # fired (CNR keeps it with heap_handle == 0)
        ages = enabled_ages(s, 0.3)
        remaining = sort([k for (k, _) in ages])
        @test remaining == [:c, :d]
        @test !any(k === :a || k === :b for (k, _) in ages)
    end

    # Pin the CNR internal-state claim directly: the fired/disabled clocks are
    # still in transition_entry (for the Anderson/Gibson-Bruck draw reuse) but
    # with heap_handle == 0, and enabled_ages must not surface them.
    rng = Xoshiro(5051)
    s = CombinedNextReaction{Symbol,Float64}()
    for k in (:a, :b, :c)
        enable!(s, k, Weibull(1.5, 2.0), 0.0, 0.0)
    end
    disable!(s, :b, 0.2)
    fire!(s, :a, 0.3)
    @test haskey(s.transition_entry, :a) && s.transition_entry[:a].heap_handle == 0
    @test haskey(s.transition_entry, :b) && s.transition_entry[:b].heap_handle == 0
    @test [k for (k, _) in enabled_ages(s, 0.3)] == [:c]
end


@testset "estimator_verbs: retained_draw on CombinedNextReaction satisfies heap_tentative == te + invlogccdf(dist, logu) to 1e-9 for a Linear-space and a Log-space family, and for a left-shifted enabling" begin
    rng = Xoshiro(6060)
    s = CombinedNextReaction{Symbol,Float64}()
    # Weibull is a Log-space family; Uniform is a Linear-space family (see
    # sampling_space). :z is left-shifted to exercise the truncation-shift term
    # in the space normalization.
    enable!(s, :wlog, Weibull(1.7, 1.0), 0.0, 0.0)
    enable!(s, :ulin, Uniform(0.0, 5.0), 0.0, 0.0)
    enable!(s, :zshift, Weibull(2.0, 1.5), -0.4, 0.0)

    checks = ((:wlog, Weibull(1.7, 1.0)),
              (:ulin, Uniform(0.0, 5.0)),
              (:zshift, Weibull(2.0, 1.5)))
    for (k, dist) in checks
        rd = retained_draw(s, k)
        tentative = s[k]                 # the heap's stored tentative firing time
        @test rd.te == (k === :zshift ? -0.4 : 0.0)
        @test 0.0 <= rd.u <= 1.0
        @test isapprox(rd.u, exp(rd.logu); atol=1e-12)
        @test isapprox(tentative, rd.te + invlogccdf(dist, rd.logu); atol=1e-9)
        # And equivalently through the total-lifetime survival identity.
        @test isapprox(rd.logu, logccdf(dist, tentative - rd.te); atol=1e-9)
    end
end


@testset "estimator_verbs: the post-force next-firing law agrees between FirstReaction and CombinedNextReaction and matches the quadrature oracle — winner probabilities, mean next-firing time, and a joint winner-by-time moment, over 20,000 replications" begin
    # Three survivor clocks including one left-shifted (te = -0.4), plus a chosen
    # clock :y forced at tstar = 0.6. After the force each survivor must be
    # distributed by its lifetime conditioned on survival past tstar, on BOTH the
    # conditional-redraw backend (which redraws at next) and the scheduling
    # backend (which keeps or redraws its stored schedule).
    survivors = [(key=:x, dist=Weibull(1.7, 1.0), te=0.0),
                 (key=:v, dist=Exponential(1.2), te=0.0),
                 (key=:z, dist=Weibull(2.0, 1.5), te=-0.4)]
    dy = Exponential(0.8)
    tstar = 0.6
    oracle_px = ev_win_probability(survivors, 1, tstar)
    oracle_pz = ev_win_probability(survivors, 3, tstar)
    oracle_mt = ev_mean_time(survivors, tstar)
    oracle_jx = ev_joint_moment(survivors, 1, tstar)
    nrep = 20_000

    stats = Dict{String,NamedTuple}()
    for (label, B) in (("FirstReaction", FirstReaction),
                       ("CombinedNextReaction", CombinedNextReaction))
        xw = Vector{Float64}(undef, nrep)
        zw = Vector{Float64}(undef, nrep)
        times = Vector{Float64}(undef, nrep)
        for r in 1:nrep
            # The sampler owns its randomness, so independent replications come
            # from distinct per-rep seeds (a fixed offset keeps the run reproducible).
            s = B{Symbol,Float64}(70_707 + r)
            for c in survivors
                enable!(s, c.key, c.dist, c.te, 0.0)
            end
            enable!(s, :y, dy, 0.0, 0.0)
            force_fire!(s, :y, tstar)
            t2, w2 = next(s, tstar)
            xw[r] = w2 === :x ? 1.0 : 0.0
            zw[r] = w2 === :z ? 1.0 : 0.0
            times[r] = t2
        end
        joint = times .* xw
        st = (px=mean(xw), sx=std(xw) / sqrt(nrep),
              pz=mean(zw), sz=std(zw) / sqrt(nrep),
              mt=mean(times), smt=std(times) / sqrt(nrep),
              mj=mean(joint), sj=std(joint) / sqrt(nrep))
        stats[label] = st
        for (est, se, oracle) in ((st.px, st.sx, oracle_px),
                                  (st.pz, st.sz, oracle_pz),
                                  (st.mt, st.smt, oracle_mt),
                                  (st.mj, st.sj, oracle_jx))
            @test abs(est - oracle) < 4 * se
            @test se < abs(oracle) / 5    # the comparison is not vacuous
        end
        # Every survivor was conditioned on surviving past tstar, so no next
        # firing may precede it — a broken repair would show up here first.
        @test all(times .> tstar)
    end
    a, b = stats["FirstReaction"], stats["CombinedNextReaction"]
    @test abs(a.px - b.px) < 4 * hypot(a.sx, b.sx)
    @test abs(a.pz - b.pz) < 4 * hypot(a.sz, b.sz)
    @test abs(a.mt - b.mt) < 4 * hypot(a.smt, b.smt)
    @test abs(a.mj - b.mj) < 4 * hypot(a.sj, b.sj)
end


@testset "estimator_verbs: on CombinedNextReaction a loser whose schedule exceeds tstar keeps its exact stored time (keep-if-later), while a loser whose schedule was passed is redrawn beyond tstar with its retained-draw identity intact (redraw-if-passed)" begin
    rng = Xoshiro(8080)
    s = CombinedNextReaction{Symbol,Float64}()
    # A spread of survivors so their scheduled times straddle a chosen tstar,
    # plus a separate clock :y to force. Peeking at the schedules to pick tstar
    # is legitimate HERE because this test certifies the mechanism (keep vs
    # redraw), not the statistical unbiasedness the precondition governs.
    survdists = Dict(:a => Weibull(1.7, 1.0), :b => Exponential(1.0),
                     :c => Weibull(2.0, 1.5), :d => Exponential(2.0),
                     :e => Uniform(0.0, 4.0))
    for (k, d) in survdists
        enable!(s, k, d, 0.0, 0.0)
    end
    enable!(s, :y, Exponential(0.5), 0.0, 0.0)

    sched_before = Dict(k => s[k] for k in keys(survdists))
    surv_sorted = sort(collect(values(sched_before)))
    tstar = surv_sorted[3]   # median-ish: some survivors below, some above
    kept = [k for (k, t) in sched_before if t > tstar]
    passed = [k for (k, t) in sched_before if t <= tstar]
    @test !isempty(kept) && !isempty(passed)   # both branches are exercised

    force_fire!(s, :y, tstar)

    for k in kept
        # keep-if-later: the stored schedule is untouched, to the bit.
        @test s[k] == sched_before[k]
    end
    for k in passed
        # redraw-if-passed: a fresh time strictly beyond tstar, and the
        # retained-draw identity still holds against the new schedule.
        @test s[k] > tstar
        rd = retained_draw(s, k)
        @test isapprox(s[k], rd.te + invlogccdf(survdists[k], rd.logu); atol=1e-9)
    end
    @test !isenabled(s, :y)
end


@testset "estimator_verbs: the capability traits report the right support per sampler, and every estimator verb throws an ArgumentError naming the missing trait when called on a sampler that does not support it" begin
    # The three aged, scheduling-or-redraw backends support enabled_ages and
    # force_fire!; only CombinedNextReaction retains a draw. The traits dispatch
    # on the bare (UnionAll) type as well as a fully parameterized one.
    for B in (FirstReaction, FirstToFire, CombinedNextReaction)
        @test supports_enabled_ages(B)
        @test supports_force(B)
    end
    @test supports_retained_draw(CombinedNextReaction)
    @test supports_retained_draw(CombinedNextReaction{Int,Float64})
    @test !supports_retained_draw(FirstReaction)
    @test !supports_retained_draw(FirstToFire)

    # The exponential-only samplers keep no enabling-time table, so age is
    # undefined and none of the estimator verbs apply.
    for B in (DirectCall, MultipleDirect, RSSA, PSSACR, Petri)
        @test !supports_enabled_ages(B)
        @test !supports_force(B)
        @test !supports_retained_draw(B)
    end
    # The trait also answers on an instance, forwarding to the type.
    @test supports_force(CombinedNextReaction{Int,Float64}())
    @test !supports_force(DirectCall{Int,Float64}())

    rng = Xoshiro(9090)
    dc = DirectCall{Symbol,Float64}()
    ftf = FirstToFire{Symbol,Float64}()
    enable!(ftf, :x, Exponential(1.0), 0.0, 0.0)
    # Unsupported verbs fail cleanly at the boundary, not with a MethodError.
    @test_throws ArgumentError enabled_ages(dc, 0.0)
    @test_throws ArgumentError force_fire!(dc, :x, 0.5)
    @test_throws ArgumentError retained_draw(dc, :x)
    @test_throws ArgumentError retained_draw(ftf, :x)   # scheduling, but no retained uniform
end


@testset "estimator_verbs: after a forced firing the sampler stays a valid live sampler — the fired clock is gone, the losers remain enabled, and normal next/fire! continues on all three aged backends" begin
    for (label, B) in EV_AGED_BACKENDS
        rng = Xoshiro(1212)
        s = B{Symbol,Float64}()
        for k in (:p, :q, :r)
            enable!(s, k, Weibull(1.5, 1.0), 0.0, 0.0)
        end
        enable!(s, :y, Exponential(0.7), 0.0, 0.0)

        force_fire!(s, :y, 0.5)
        @test !isenabled(s, :y)
        @test all(isenabled(s, k) for k in (:p, :q, :r))
        @test sort([k for (k, _) in enabled_ages(s, 0.5)]) == [:p, :q, :r]

        # The sampler keeps working: next returns one of the survivors past the
        # forced time, and firing it leaves the rest enabled.
        t2, w2 = next(s, 0.5)
        @test w2 in (:p, :q, :r)
        @test t2 > 0.5
        fire!(s, w2, t2)
        @test !isenabled(s, w2)
        @test length([k for (k, _) in enabled_ages(s, t2)]) == 2
    end
end


@testset "estimator_verbs: a context-level force_fire! advances time and fans out to watchers, so an attached TrajectoryRecorder records the forced firing with the back-calculated retained draw u == ccdf(dist, tstar - te)" begin
    using CompetingClocks: SamplingContext, SamplerBuilder, NextReactionMethod,
        with_recorder, close_record!, recorded_firings, force_fire!

    dist_y = Exponential(0.8)
    ctx = SamplingContext(SamplerBuilder(Symbol, Float64; method=NextReactionMethod()),
                          Xoshiro(31_337))
    ctx, rec = with_recorder(ctx)
    enable!(ctx, :x, Weibull(1.7, 1.0))
    enable!(ctx, :y, dist_y)
    enable!(ctx, :z, Weibull(2.0, 1.5))

    tstar = 0.55
    force_fire!(ctx, :y, tstar)
    @test time(ctx) == tstar

    close_record!(rec, 10.0)
    firings = recorded_firings(rec)
    @test length(firings) == 1
    fr = firings[1]
    @test fr.clock === :y
    @test fr.when == tstar
    # For a forced firing the recorder's back-calculated uniform is still the
    # total-lifetime survival ccdf(dist, tstar - te), so a replay cannot tell an
    # imposed firing from a raced one.
    @test isapprox(fr.u, ccdf(dist_y, tstar - fr.te); atol=1e-12)
    @test isapprox(fr.when, fr.te + invlogccdf(dist_y, fr.logu); atol=1e-9)
    # The losers remain enabled and the sampler advanced past the force.
    @test isenabled(ctx, :x) && isenabled(ctx, :z) && !isenabled(ctx, :y)
end
