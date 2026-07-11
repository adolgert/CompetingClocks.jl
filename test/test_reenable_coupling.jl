# Tests for the re-evaluation verb and the sampler-owned coupling (Milestone 5,
# guarantee G6, revised): `reenable!(sampler, key, dist, te, when)` re-evaluates
# a still-enabled clock's distribution, and HOW the change is realized pathwise
# — :carry (deterministic, IPA-safe) or :redraw (redraw-at-change) — is a
# construction-time property of the sampler (`coupling=` keyword, read back with
# `CompetingClocks.coupling`), not a per-call argument. Both couplings agree in
# law; they differ only pathwise. Plain @testset (not @safetestset) so the
# file-scope oracles below are shared.
#   julia --project=test test/runtests.jl "reenable"
using Test
using CompetingClocks
using CompetingClocks: FirstReaction, FirstToFire, CombinedNextReaction,
    DirectCall, MultipleDirect, RSSA, PSSACR, Petri, SamplerChoice,
    supports_carry, validate_coupling, coupling, reenable!,
    enable!, disable!, fire!, next, isenabled,
    clone, similar_sampler, copy_clocks!,
    SamplingContext, SamplerBuilder, NextReactionMethod, FirstToFireMethod,
    with_recorder, close_record!, recorded_firings, TrajectoryWatcher,
    pathloglikelihood
using Distributions: Exponential, Weibull, Gamma, Uniform, rate, scale,
    logccdf, logpdf, invlogccdf, ccdf, cdf
using Random: Xoshiro
using Statistics: mean, std


# --- analytic oracle: single clock, piecewise-exponential hazard -------------
# A clock enabled at te=0 with rate λ0, whose rate switches to λ1 at absolute
# time τ (age τ, since te=0). Firing-time survival is S(t)=exp(-λ0 t) for t<τ and
# exp(-λ0 τ - λ1 (t-τ)) for t≥τ. Its mean is the tail integral of S.
piecewise_exp_mean(λ0, λ1, τ) = (1 - exp(-λ0 * τ)) / λ0 + exp(-λ0 * τ) / λ1


# MultipleDirect needs a chooser at construction; one bucket is enough to
# exercise its constructor-time coupling validation.
struct ReenableOneBucket <: SamplerChoice{Symbol,Symbol} end


# Drive one single-clock piecewise-exponential firing on a bare sampler: enable
# at rate λ0, and if the (old-law) draw survives to τ, re-evaluate to rate λ1 at
# age τ, then read the firing time. The coupling is baked into the sampler at
# construction. Returns the firing time.
function drive_single(B, λ0, λ1, τ, cpl, seed)
    s = B{Symbol,Float64}(seed; coupling=cpl)
    enable!(s, :A, Exponential(1 / λ0), 0.0, 0.0)
    t, _ = next(s, 0.0)
    t < τ && return t                      # fired before the change, old law
    reenable!(s, :A, Exponential(1 / λ1), 0.0, τ)  # te=0 keeps the age
    t2, _ = next(s, τ)
    return t2
end


# Drive one two-clock race: A is piecewise-exponential (rate λA0 then λA1 at τ),
# B is constant-rate λB. If nothing fires before τ, re-evaluate A at age τ and
# continue. Returns (winner_is_A, firing_time).
function drive_race(B, λA0, λA1, λB, τ, cpl, seed)
    s = B{Symbol,Float64}(seed; coupling=cpl)
    enable!(s, :A, Exponential(1 / λA0), 0.0, 0.0)
    enable!(s, :B, Exponential(1 / λB), 0.0, 0.0)
    t, w = next(s, 0.0)
    t < τ && return (w === :A, t)
    reenable!(s, :A, Exponential(1 / λA1), 0.0, τ)
    t2, w2 = next(s, τ)
    return (w2 === :A, t2)
end


@testset "reenable: the supports_carry trait is true exactly for the backends that retain an in-flight draw to rebase (CombinedNextReaction, FirstToFire) and false for the redraw-only and memoryless backends" begin
    # Carry needs a retained draw to map through the change. CNR keeps the
    # survival uniform; FTF reconstructs it from the stored firing time. The
    # exponential-only samplers keep no draw, and FirstReaction redraws at every
    # next, so none of them can carry.
    @test supports_carry(CombinedNextReaction)
    @test supports_carry(FirstToFire)
    @test supports_carry(CombinedNextReaction{Int,Float64})
    @test !supports_carry(FirstReaction)
    for Bexp in (DirectCall, MultipleDirect, RSSA, PSSACR, Petri)
        @test !supports_carry(Bexp)
    end
    # The trait answers on an instance, forwarding to the type.
    @test supports_carry(FirstToFire{Int,Float64}())
    @test !supports_carry(DirectCall{Int,Float64}())
end


@testset "reenable coupling: constructing any sampler with an unknown coupling symbol throws, constructing a redraw-only sampler with coupling=:carry throws naming supports_carry, and the carry-capable samplers accept both couplings" begin
    # The constructor is the single validation point: a bad coupling must fail
    # at construction, not at the first reenable!.
    @test_throws ArgumentError CombinedNextReaction{Symbol,Float64}(coupling=:bogus)
    @test_throws ArgumentError FirstToFire{Symbol,Float64}(coupling=:bogus)
    @test_throws ArgumentError FirstReaction{Symbol,Float64}(coupling=:bogus)
    @test_throws ArgumentError DirectCall{Symbol,Float64}(coupling=:bogus)

    # The redraw-only samplers accept :redraw (the behavior they already have)
    # and reject :carry with an error that names the missing capability.
    builders = (
        ("FirstReaction", (; kw...) -> FirstReaction{Symbol,Float64}(; kw...)),
        ("DirectCall", (; kw...) -> DirectCall{Symbol,Float64}(; kw...)),
        ("MultipleDirect", (; kw...) -> MultipleDirect{Symbol,Symbol,Float64}(ReenableOneBucket(); kw...)),
        ("RSSA", (; kw...) -> RSSA{Symbol,Float64}(; kw...)),
        ("PSSACR", (; kw...) -> PSSACR{Symbol,Float64}(; kw...)),
        ("Petri", (; kw...) -> Petri{Symbol,Float64}(1.0; kw...)),
    )
    for (name, build) in builders
        # :redraw builds fine and is what the accessor reports.
        s = build(coupling=:redraw)
        @test coupling(s) === :redraw
        # :carry fails at construction with the trait named in the message.
        err = try
            build(coupling=:carry)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("supports_carry", err.msg)
    end

    # The carry-capable samplers accept both couplings.
    @test coupling(CombinedNextReaction{Symbol,Float64}(coupling=:carry)) === :carry
    @test coupling(CombinedNextReaction{Symbol,Float64}(coupling=:redraw)) === :redraw
    @test coupling(FirstToFire{Symbol,Float64}(coupling=:carry)) === :carry
    @test coupling(FirstToFire{Symbol,Float64}(coupling=:redraw)) === :redraw
end


@testset "reenable coupling: the coupling accessor defaults to :carry on the carry-capable samplers, answers :redraw generically for memoryless samplers, and a SamplingContext forwards to its inner sampler" begin
    # Default construction restores the historical silent-carry behavior of the
    # scheduling backends as the explicit default.
    @test coupling(CombinedNextReaction{Symbol,Float64}()) === :carry
    @test coupling(FirstToFire{Symbol,Float64}()) === :carry
    # The generic SSA fallback serves every sampler without a field.
    @test coupling(DirectCall{Symbol,Float64}()) === :redraw
    @test coupling(FirstReaction{Symbol,Float64}()) === :redraw
    @test coupling(Petri{Symbol,Float64}()) === :redraw

    ctx_carry = SamplingContext(SamplerBuilder(Symbol, Float64;
            method=NextReactionMethod()), Xoshiro(5))
    @test coupling(ctx_carry) === :carry
    ctx_redraw = SamplingContext(SamplerBuilder(Symbol, Float64;
            method=NextReactionMethod(coupling=:redraw)), Xoshiro(5))
    @test coupling(ctx_redraw) === :redraw
end


@testset "reenable coupling: the NextReactionMethod and FirstToFireMethod specs carry a coupling keyword, validate it at spec construction, and build samplers whose coupling field matches" begin
    # The spec is how framework users select samplers, so it must forward the
    # coupling to the sampler constructor and fail as early as the constructor
    # would.
    @test coupling(NextReactionMethod()(Symbol, Float64)) === :carry
    @test coupling(NextReactionMethod(coupling=:redraw)(Symbol, Float64)) === :redraw
    @test coupling(FirstToFireMethod()(Symbol, Float64)) === :carry
    @test coupling(FirstToFireMethod(coupling=:redraw)(Symbol, Float64)) === :redraw
    @test_throws ArgumentError NextReactionMethod(coupling=:bogus)
    @test_throws ArgumentError FirstToFireMethod(coupling=:bogus)
end


@testset "reenable coupling: clone, similar_sampler, and copy_clocks! all preserve the coupling field on CombinedNextReaction and FirstToFire" begin
    # The coupling is part of the sampler's identity: a clone that silently
    # flipped back to the default would decouple a pathwise experiment.
    for B in (CombinedNextReaction, FirstToFire)
        s = B{Symbol,Float64}(31; coupling=:redraw)
        enable!(s, :A, Weibull(1.5, 1.0), 0.0, 0.0)
        @test coupling(clone(s)) === :redraw
        @test coupling(similar_sampler(s)) === :redraw
        dst = B{Symbol,Float64}(32)             # built :carry by default
        copy_clocks!(dst, s)
        @test coupling(dst) === :redraw
    end
end


@testset "reenable: an unknown clock is rejected, and on a memoryless sampler reenable! is the redraw re-evaluation that replaces the stored rate" begin
    ftf = FirstToFire{Symbol,Float64}()
    enable!(ftf, :x, Exponential(1.0), 0.0, 0.0)
    # Not-enabled clock.
    @test_throws ArgumentError reenable!(ftf, :missing, Exponential(1.0), 0.0, 0.5)

    # On an exponential-only backend the generic reenable! just replaces the
    # memoryless rate — redraw is the only behavior such a sampler has.
    dc = DirectCall{Symbol,Float64}()
    enable!(dc, :x, Exponential(1.0), 0.0, 0.0)
    reenable!(dc, :x, Exponential(0.25), 0.0, 0.5)
    @test dc[:x] == rate(Exponential(0.25))
end


@testset "reenable: a :carry sampler re-evaluated with an unchanged distribution leaves the schedule bit-for-bit unchanged while a :redraw sampler generally moves it, on both CombinedNextReaction and FirstToFire" begin
    # This is the structural signature of the two couplings: carry is a
    # deterministic map that is the identity when the distribution does not
    # change; redraw consumes a fresh number from the clock's stream, so even an
    # identical distribution generally reschedules.
    for B in (CombinedNextReaction, FirstToFire)
        dist = Weibull(1.6, 1.3)
        s = B{Symbol,Float64}(2024; coupling=:carry)
        enable!(s, :A, dist, 0.0, 0.0)
        sched0 = s[:A]
        reenable!(s, :A, dist, 0.0, 0.4)
        @test s[:A] == sched0                       # carry + same dist = identity

        s2 = B{Symbol,Float64}(2024; coupling=:redraw)
        enable!(s2, :A, dist, 0.0, 0.0)
        sched0b = s2[:A]
        reenable!(s2, :A, dist, 0.0, 0.4)
        @test s2[:A] != sched0b                      # redraw draws fresh
        @test s2[:A] > 0.4                           # and past the change time
    end
end


@testset "reenable: under a :carry sampler the new schedule is a deterministic, continuous function of a small rate change (moves O(epsilon)), while a :redraw sampler jumps, on CombinedNextReaction and FirstToFire" begin
    # Pathwise smoke check that pins the difference the milestone exists for. A
    # tiny change in the new rate moves the carried schedule by O(epsilon)
    # (continuous, RNG-free), whereas redraw's schedule is a fresh draw unrelated
    # to epsilon.
    for B in (CombinedNextReaction, FirstToFire)
        τ = 0.3
        function carried_sched(newrate)
            s = B{Symbol,Float64}(777; coupling=:carry)
            enable!(s, :A, Exponential(1 / 0.8), 0.0, 0.0)  # rate 0.8
            reenable!(s, :A, Exponential(1 / newrate), 0.0, τ)
            return s[:A]
        end
        base = carried_sched(1.5)
        d1 = abs(carried_sched(1.5 + 1e-4) - base)
        d2 = abs(carried_sched(1.5 + 2e-4) - base)
        @test d1 < 1e-2                              # small change → small move
        @test isapprox(d2 / d1, 2.0; rtol=1e-2)      # linear in epsilon (O(ε))

        # Redraw does not vary continuously with the new rate at fixed seed: it is
        # a fresh draw from the clock's stream. Confirm it differs from carry.
        s = B{Symbol,Float64}(777; coupling=:redraw)
        enable!(s, :A, Exponential(1 / 0.8), 0.0, 0.0)
        reenable!(s, :A, Exponential(1 / 1.5), 0.0, τ)
        @test s[:A] != base
    end
end


@testset "reenable: for exponential old and new distributions, a :carry sampler reduces EXACTLY to Gibson-Bruck rate-ratio scaling of the remaining time to relative tolerance 1e-12, on CombinedNextReaction and FirstToFire" begin
    # The pin the design doc names: carry of an exponential clock scales its
    # remaining time by rate_old/rate_new. Anything else means the carry map is
    # not the memoryless-consistent one.
    for B in (CombinedNextReaction, FirstToFire)
        rate_old, rate_new = 0.7, 2.3
        τ = 0.15
        s = B{Symbol,Float64}(4242; coupling=:carry)
        enable!(s, :A, Exponential(1 / rate_old), 0.0, 0.0)
        sched_old = s[:A]
        @test sched_old > τ                          # the clock is still pending at τ
        reenable!(s, :A, Exponential(1 / rate_new), 0.0, τ)
        sched_new = s[:A]
        @test isapprox(sched_new - τ, (rate_old / rate_new) * (sched_old - τ); rtol=1e-12)
    end
end


@testset "reenable: FirstReaction treats reenable! as the moot redraw-everything re-evaluation, replacing the stored distribution" begin
    # FirstReaction retains no schedule and redraws at every next, so its only
    # possible re-evaluation is to replace the distribution; no coupling field
    # is stored.
    s = FirstReaction{Symbol,Float64}(99)
    enable!(s, :A, Exponential(2.0), 0.0, 0.0)
    reenable!(s, :A, Gamma(2.0, 0.5), 0.0, 0.5)
    @test isenabled(s, :A)
    @test s.enabled[:A].distribution == Gamma(2.0, 0.5)
end


@testset "reenable: on a single clock with a piecewise-exponential hazard, samplers built with either coupling on both scheduling backends and FirstReaction reproduce the analytic piecewise firing-time mean within 4 standard errors, over 40,000 replications" begin
    λ0, λ1, τ = 0.6, 1.8, 0.7
    oracle = piecewise_exp_mean(λ0, λ1, τ)
    nrep = 40_000
    cases = (
        ("CNR/carry", CombinedNextReaction, :carry),
        ("CNR/redraw", CombinedNextReaction, :redraw),
        ("FTF/carry", FirstToFire, :carry),
        ("FTF/redraw", FirstToFire, :redraw),
        ("FR/redraw", FirstReaction, :redraw),
    )
    for (label, B, cpl) in cases
        times = Vector{Float64}(undef, nrep)
        for r in 1:nrep
            times[r] = drive_single(B, λ0, λ1, τ, cpl, 100_000 + r)
        end
        se = std(times) / sqrt(nrep)
        @test abs(mean(times) - oracle) < 4 * se     # matches the piecewise law
        @test se < oracle / 20                        # the comparison is not vacuous
    end
end


@testset "reenable: in a two-clock race with a deterministic mid-flight re-evaluation of one clock, a :carry and a :redraw CombinedNextReaction agree with each other and with FirstReaction on the winner probability and mean firing time, within 4 pooled standard errors, over 40,000 replications" begin
    λA0, λA1, λB, τ = 0.5, 2.0, 1.0, 0.6
    nrep = 40_000
    function run_case(B, cpl, base)
        winA = Vector{Float64}(undef, nrep)
        times = Vector{Float64}(undef, nrep)
        for r in 1:nrep
            a, t = drive_race(B, λA0, λA1, λB, τ, cpl, base + r)
            winA[r] = a ? 1.0 : 0.0
            times[r] = t
        end
        return (pA=mean(winA), spA=std(winA) / sqrt(nrep),
                mt=mean(times), smt=std(times) / sqrt(nrep))
    end
    carry = run_case(CombinedNextReaction, :carry, 200_000)
    redraw = run_case(CombinedNextReaction, :redraw, 300_000)
    fr = run_case(FirstReaction, :redraw, 400_000)

    for (a, b) in ((carry, redraw), (carry, fr), (redraw, fr))
        @test abs(a.pA - b.pA) < 4 * hypot(a.spA, b.spA)
        @test abs(a.mt - b.mt) < 4 * hypot(a.smt, b.smt)
    end
    # The comparisons must have teeth: the winner probability is well inside (0,1)
    # and its standard error is small.
    @test 0.1 < carry.pA < 0.9
    @test carry.spA < 0.02
end


@testset "reenable: a TrajectoryWatcher forward-accumulates exactly the hand-computed two-segment piecewise log-likelihood when a clock's distribution changes mid-flight" begin
    # The likelihood layer must see a re-enable as a segment boundary: close the
    # OLD segment's survival, open the NEW one. The watcher's inherited enable!
    # disables-then-restores, and its disable! banks the old-segment survival, so
    # this happens with no special reenable! code on the watcher itself.
    dist_old = Weibull(1.4, 1.0)
    dist_new = Gamma(2.0, 0.7)
    τ = 0.9
    tfire = 2.1

    tw = TrajectoryWatcher{Symbol,Float64}()
    enable!(tw, :A, dist_old, 0.0, 0.0)
    enable!(tw, :A, dist_new, 0.0, τ)     # the mid-flight change (te=0 keeps age)
    fire!(tw, :A, tfire)

    # Hand piecewise value: survive the old law over [0, τ], then fire under the
    # new law with its own survival back to τ subtracted.
    hand = logccdf(dist_old, τ) + logpdf(dist_new, tfire) - logccdf(dist_new, τ)
    @test isapprox(tw.loglikelihood, hand; atol=1e-12)
end


@testset "reenable: a context-level reenable! keeping the clock's age routes the segment boundary to the path-likelihood watcher, so the accumulated path log-likelihood equals the hand-computed piecewise value" begin
    # A trigger clock advances context time to τ; then we re-evaluate the tracked
    # clock keeping its age (relative shift = original_te - now = -τ) and fire it.
    dist_old = Weibull(1.5, 1.2)
    dist_new = Exponential(0.8)
    τ = 0.5
    tfire = 1.7

    ctx = SamplingContext(SamplerBuilder(Symbol, Float64;
            method=NextReactionMethod(), path_likelihood=true), Xoshiro(11))
    enable!(ctx, :A, dist_old)            # te = 0
    enable!(ctx, :trigger, Exponential(1.0))
    fire!(ctx, :trigger, τ)               # advances ctx.time to τ
    reenable!(ctx, :A, dist_new, -τ)      # keep age: te stays 0; sampler carries
    fire!(ctx, :A, tfire)

    # After :A fires, the only remaining survival is none; path loglik at tfire is
    # the closed 2-segment value for :A (trigger contributed its own fire term).
    ll = pathloglikelihood(ctx, tfire)
    seg_A = logccdf(dist_old, τ) + logpdf(dist_new, tfire) - logccdf(dist_new, τ)
    trigger_term = logpdf(Exponential(1.0), τ)
    @test isapprox(ll, seg_A + trigger_term; atol=1e-10)
end


@testset "reenable: a context-level reenable! updates the TrajectoryRecorder's live table so the fired clock's back-calculated uniform uses the CURRENT segment's distribution, not the retired one" begin
    # The recorder back-calculates u == ccdf(dist, when - te) from its live
    # (dist, te) table. A mid-flight change must update that table, or the u would
    # be computed against the retired distribution. The sampler is built with
    # coupling=:redraw to exercise the redraw path through the context.
    dist_old = Weibull(1.3, 1.0)
    dist_new = Exponential(0.6)
    τ = 0.4
    tfire = 1.3

    ctx = SamplingContext(SamplerBuilder(Symbol, Float64;
            method=NextReactionMethod(coupling=:redraw)), Xoshiro(22))
    ctx, rec = with_recorder(ctx)
    enable!(ctx, :A, dist_old)
    enable!(ctx, :trigger, Exponential(1.0))
    fire!(ctx, :trigger, τ)
    reenable!(ctx, :A, dist_new, -τ)      # keep age; the sampler redraws
    # :A fires at tfire (drive it directly; the exact time is what we record).
    fire!(ctx, :A, tfire)

    close_record!(rec, 5.0)
    firings = recorded_firings(rec)
    frA = only(f for f in firings if f.clock === :A)
    @test frA.distribution == dist_new                         # current segment
    @test frA.te == 0.0                                        # age kept
    @test isapprox(frA.u, ccdf(dist_new, tfire - 0.0); atol=1e-12)
    @test isapprox(frA.when, frA.te + invlogccdf(dist_new, frA.logu); atol=1e-9)
end
