using SafeTestsets


@safetestset CombinedNextReactionSmoke = "CombinedNextReaction reaction does basic things" begin
    using Distributions
    using Random
    using CompetingClocks: CombinedNextReaction, next, enable!, disable!, reset!

    rng = MersenneTwister(349827)
    for i in 1:100
        sampler = CombinedNextReaction{String,Float64}()
        @test next(sampler, 3.0)[2] === nothing
        enable!(sampler, "walk home", Exponential(1.5), 0.0, 0.0)
        @test next(sampler, 3.0)[2] == "walk home"
        enable!(sampler, "run", Gamma(1, 3), 0.0, 0.0)
        @test next(sampler, 3.0)[2] ∈ ["walk home", "run"]
        enable!(sampler, "walk to sandwich shop", Weibull(2, 1), 0.0, 0.0)
        @test next(sampler, 3.0)[2] ∈ ["walk home", "run", "walk to sandwich shop"]
        disable!(sampler, "walk to sandwich shop", 1.7)
        @test next(sampler, 3.0)[2] ∈ ["walk home", "run"]
        reset!(sampler)
    end
end

@safetestset CombinedNextReaction_interface = "CombinedNextReaction basic interface" begin
    using CompetingClocks
    using CompetingClocks: CombinedNextReaction, jitter!, sampling_space
    using Distributions
    using Random: Xoshiro

    sampler = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[1]
    @test keytype(sampler) <: Int64

    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[1] == 7.9

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "1")

    disable!(sampler, 1, 0.0)

    # A disabled clock keeps a retained-survival entry (heap_handle == 0) but
    # has no schedule; getindex treats it as absent, matching haskey.
    @test_throws KeyError sampler[1]
    @test sampler[2] == 12.3

end


@safetestset CombinedNextReaction_pair_keys = "CombinedNextReaction pair keys" begin
    using CompetingClocks
    using CompetingClocks: CombinedNextReaction, jitter!, sampling_space
    using Distributions
    using Random: Xoshiro

    # Make a key that looks like those in Gen.jl.
    KeyType = Pair{Symbol,Int64}
    sampler = CombinedNextReaction{KeyType,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[:S=>1]
    @test keytype(sampler) <: KeyType

    for (clock, when_fire) in [
        (:S => 1, 7.9), (:S => 2, 12.3), (:I => 3, 3.7), (:I => 4, 0.00013), (:S => 5, 0.2)
    ]
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[:S=>1] == 7.9

    @test haskey(sampler, :S => 1)
    @test !haskey(sampler, :I => 17)

    disable!(sampler, :S => 1, 0.0)

    # Same absence semantics as the Int-keyed test above: KeyError, not a
    # BoundsError from indexing the heap at the retained entry's zero handle.
    @test_throws KeyError sampler[:S=>1]
    @test sampler[:S=>2] == 12.3

end

@safetestset CombinedNextReaction_copy = "CombinedNextReaction copy" begin
    using CompetingClocks
    using CompetingClocks: CombinedNextReaction, jitter!, sampling_space
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    dst = clone(src)
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0)
    enable!(src, 38, Exponential(), 0.0, 0.0)
    enable!(dst, 29, Exponential(), 0.0, 0.0)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(src) == 2
    @test length(dst) == 2
    # Changing src doesn't change dst.
    enable!(src, 48, Exponential(), 0.0, 0.0)
    @test length(src) == 3
    @test length(dst) == 2
    # Changing dst doesn't change src.
    enable!(dst, 49, Exponential(), 0.0, 0.0)
    @test length(src) == 3
    @test length(dst) == 3
end


@safetestset CombinedNextReaction_set = "CombinedNextReaction set" begin
    using CompetingClocks
    using CompetingClocks: CombinedNextReaction, jitter!, sampling_space
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0)
    enable!(src, 38, Exponential(), 0.0, 0.0)
    enable!(src, 29, Exponential(), 0.0, 0.0)
    disable!(src, 37, 0.1)

    enabled_set = enabled(src)
    @test 38 in enabled_set
    @test 37 ∉ enabled_set
    @test length(enabled_set) == 2
    @test eltype(enabled_set) == Int64
    b = Set([x for x in enabled_set])
    @assert b == Set([38, 29])
end


@safetestset CombinedNextReaction_sampling_space = "CombinedNextReaction sampling_space" begin
    using CompetingClocks: sampling_space, LinearSampling, LogSampling
    using Distributions

    # Test that sampling_space returns correct types for various distributions
    # LinearSampling distributions
    @test sampling_space(Arcsine) == LinearSampling
    @test sampling_space(Beta) == LinearSampling
    @test sampling_space(BetaPrime) == LinearSampling
    @test sampling_space(Biweight) == LinearSampling
    @test sampling_space(Cauchy) == LinearSampling
    @test sampling_space(Cosine) == LinearSampling
    @test sampling_space(Epanechnikov) == LinearSampling
    @test sampling_space(Frechet) == LinearSampling
    @test sampling_space(GeneralizedPareto) == LinearSampling
    @test sampling_space(Gumbel) == LinearSampling
    @test sampling_space(InverseGamma) == LinearSampling
    @test sampling_space(InverseGaussian) == LinearSampling
    @test sampling_space(Kolmogorov) == LinearSampling
    @test sampling_space(Kumaraswamy) == LinearSampling
    @test sampling_space(Levy) == LinearSampling
    @test sampling_space(Logistic) == LinearSampling
    @test sampling_space(LogitNormal) == LinearSampling
    @test sampling_space(LogNormal) == LinearSampling
    @test sampling_space(Normal) == LinearSampling
    @test sampling_space(NormalCanon) == LinearSampling
    @test sampling_space(Pareto) == LinearSampling
    @test sampling_space(PGeneralizedGaussian) == LinearSampling
    @test sampling_space(Rayleigh) == LinearSampling
    @test sampling_space(SymTriangularDist) == LinearSampling
    @test sampling_space(Triweight) == LinearSampling
    @test sampling_space(Uniform) == LinearSampling

    # LogSampling distributions
    @test sampling_space(Erlang) == LogSampling
    @test sampling_space(Exponential) == LogSampling
    @test sampling_space(Gamma) == LogSampling
    @test sampling_space(Laplace) == LogSampling
    @test sampling_space(Weibull) == LogSampling

    # Test with instances - works because Type{<:Distribution} matches parametric types
    @test sampling_space(Exponential(1.0)) == LogSampling
    @test sampling_space(Normal(0, 1)) == LinearSampling
    @test sampling_space(Gamma(2.0, 1.0)) == LogSampling
    @test sampling_space(Uniform(0, 1)) == LinearSampling
end


@safetestset CombinedNextReaction_survival_helpers = "CombinedNextReaction survival helpers" begin
    using CompetingClocks: get_survival_zero, draw_space, survival_space, invert_space
    using CompetingClocks: LinearSampling, LogSampling
    using Random: Xoshiro
    using Distributions

    # Test get_survival_zero
    @test get_survival_zero(LinearSampling) == 0.0
    @test get_survival_zero(LogSampling) == -Inf
    @test get_survival_zero(Exponential) == -Inf  # LogSampling
    @test get_survival_zero(Normal) == 0.0  # LinearSampling
    @test get_survival_zero(Exponential(1.0)) == -Inf  # Instance test

    # Test draw_space
    rng = Xoshiro(12345)
    linear_draw = draw_space(LinearSampling, rng)
    @test 0.0 <= linear_draw <= 1.0  # Uniform draw

    rng = Xoshiro(12345)
    log_draw = draw_space(LogSampling, rng)
    @test log_draw >= 0.0  # Exponential draw

    # Test survival_space
    dist = Exponential(1.0)
    @test survival_space(Exponential, dist, 0.5) == logccdf(dist, 0.5)

    dist_normal = Normal(0, 1)
    @test survival_space(Normal, dist_normal, 0.5) == ccdf(dist_normal, 0.5)

    # Test invert_space
    dist = Exponential(1.0)
    survival = logccdf(dist, 0.5)
    @test invert_space(Exponential, dist, survival) ≈ 0.5 atol=1e-10

    dist_normal = Normal(0, 1)
    survival_linear = ccdf(dist_normal, 0.5)
    @test invert_space(Normal, dist_normal, survival_linear) ≈ 0.5 atol=1e-10
end


@safetestset CombinedNextReaction_jitter = "CombinedNextReaction jitter!" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, jitter!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(234567)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable some clocks
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0)

    # Get current next event
    t1, k1 = next(sampler, 0.0)

    # Jitter should resample all clocks
    jitter!(sampler, 0.5)

    # Get new next event - times should be different after jitter
    t2, k2 = next(sampler, 0.5)
    @test t2 >= 0.5  # Times are after the jitter point
end


@safetestset CombinedNextReaction_truncated = "CombinedNextReaction truncated sampling" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, disable!
    using Random: Xoshiro
    using Distributions: Exponential, Gamma, Weibull

    rng = Xoshiro(345678)

    # Test truncated sampling: enable with te < when
    # This exercises the truncated distribution branch in sample_shifted
    # The distribution starts at te=0 but current time when=0.5
    # So the sampler uses truncated distribution conditioned on survival past 0.5
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.5)
    t, k = next(sampler, 0.5)
    @test t >= 0.5  # Fire time must be >= current time (truncated dist guarantees this)

    # Test with LogSampling distribution (Gamma)
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.5)
    t2, k2 = next(sampler2, 0.5)
    @test t2 >= 0.5

    # Test with Weibull (also LogSampling)
    sampler3 = CombinedNextReaction{Int,Float64}()
    enable!(sampler3, 1, Weibull(2.0, 1.0), 0.0, 0.5)
    t3, k3 = next(sampler3, 0.5)
    @test t3 >= 0.5
end


@safetestset CombinedNextReaction_reenable = "CombinedNextReaction re-enable after fire" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, fire!, disable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable and fire a clock
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    t, k = next(sampler, 0.0)
    fire!(sampler, 1, t)

    # Re-enable the same clock (should exercise heap_handle > 0 but survival == 0 branch)
    enable!(sampler, 1, Exponential(1.0), t, t)
    t2, k2 = next(sampler, t)
    @test k2 == 1
    @test t2 >= t
end


@safetestset CombinedNextReaction_update_enabled = "CombinedNextReaction update enabled clock" begin
    using CompetingClocks: CombinedNextReaction, enable!, next
    using Random: Xoshiro
    using Distributions: Exponential, Gamma

    rng = Xoshiro(567890)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable a clock
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    t1, k1 = next(sampler, 0.0)

    # Update the same clock with different distribution (exercises consume_survival).
    # The re-enable time must not step past the pending firing t1: next()'s time
    # argument may never advance past a pending firing without that event being
    # fired. Survival mark moved from next! to fire! (2026-07-04): previously
    # next() zeroed the returned clock's survival, so a re-enable that stepped
    # over the pending firing (the old when=0.5) silently redrew instead of
    # reusing; with next() now pure, reuse is faithful and an out-of-contract
    # re-enable time would consume survival past the firing, so we re-enable at a
    # time strictly before t1.
    when_re = t1 / 2
    enable!(sampler, 1, Exponential(2.0), 0.0, when_re)
    t2, k2 = next(sampler, when_re)
    @test t2 >= when_re

    # Update with a LogSampling distribution to exercise that branch, again
    # re-enabling in-contract (before the pending firing t3).
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.0)
    t3, k3 = next(sampler2, 0.0)
    when_re2 = t3 / 2
    enable!(sampler2, 1, Gamma(3.0, 1.0), 0.0, when_re2)
    t4, k4 = next(sampler2, when_re2)
    @test t4 >= when_re2
end


@safetestset CombinedNextReaction_isenabled = "CombinedNextReaction isenabled" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, isenabled
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(678901)
    sampler = CombinedNextReaction{Int,Float64}()

    @test !isenabled(sampler, 1)

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    @test isenabled(sampler, 1)
    @test !isenabled(sampler, 2)

    disable!(sampler, 1, 0.5)
    @test !isenabled(sampler, 1)
end


@safetestset CombinedNextReaction_linear_distributions = "CombinedNextReaction with LinearSampling distributions" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, disable!
    using Random: Xoshiro
    using Distributions: Normal, Uniform, LogNormal, Pareto, truncated

    rng = Xoshiro(789012)

    # Test with Normal distribution (LinearSampling)
    sampler1 = CombinedNextReaction{Int,Float64}()
    enable!(sampler1, 1, truncated(Normal(5.0, 1.0), 0, Inf), 0.0, 0.0)
    t1, k1 = next(sampler1, 0.0)
    @test t1 > 0

    # Test with Uniform distribution
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Uniform(1.0, 5.0), 0.0, 0.0)
    t2, k2 = next(sampler2, 0.0)
    @test 1.0 <= t2 <= 5.0

    # Test with LogNormal
    sampler3 = CombinedNextReaction{Int,Float64}()
    enable!(sampler3, 1, LogNormal(0.0, 1.0), 0.0, 0.0)
    t3, k3 = next(sampler3, 0.0)
    @test t3 > 0

    # Test with Pareto
    sampler4 = CombinedNextReaction{Int,Float64}()
    enable!(sampler4, 1, Pareto(1.0), 0.0, 0.0)
    t4, k4 = next(sampler4, 0.0)
    @test t4 >= 1.0
end


@safetestset CombinedNextReaction_future_enable = "CombinedNextReaction enable in future" begin
    using CompetingClocks: CombinedNextReaction, enable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(890123)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable with te > when (future enabling time)
    # This exercises the else branch in sample_shifted where te >= when
    enable!(sampler, 1, Exponential(1.0), 2.0, 0.0)
    t, k = next(sampler, 0.0)
    @test t >= 2.0  # Fire time must be after the enabling time
end


@safetestset CombinedNextReaction_disable_reenable = "CombinedNextReaction disable then re-enable" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, next
    using Random: Xoshiro
    using Distributions: Exponential, Gamma

    rng = Xoshiro(901234)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable, disable, then re-enable (should use remaining survival)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)
    disable!(sampler, 1, 0.3)

    # Re-enable the disabled clock (heap_handle == 0 but has remaining survival)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.5)
    t, k = next(sampler, 0.5)
    @test k == 1
    @test t >= 0.5

    # Same test with LogSampling distribution
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.0)
    disable!(sampler2, 1, 0.3)
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.5)
    t2, k2 = next(sampler2, 0.5)
    @test k2 == 1
    @test t2 >= 0.5
end


# ======================================================================
# Cache-protection battery for CombinedNextReaction (CNR).
#
# These tests pin the Anderson/Gibson-Bruck draw-reuse machinery after the
# behavior change of 2026-07-04: the "survival mark" that committed a clock to
# firing moved from next!() to fire!(). next() is now a pure read of the firing
# queue; fire!() is what consumes a draw completely. The tests below verify:
#   (a) next() purity + zero RNG consumption
#   (b) fired path -> re-enable takes the fresh-redraw branch
#   (c) LOSER path (the fixed bug): a clock returned by next() but NOT fired
#       must, on re-enable, REUSE its draw (no RNG, analytic match)
#   (d) disable-without-fire -> re-enable reuses remaining survival (no RNG)
#   (e) consume_survival analytic checks in LinearSampling and LogSampling
#   (f) same-te/same-dist re-enable no-op (heap entry + survival unchanged, no RNG)
#   (g) MultiSampler routing regression (fire! reaches sub-sampler) + DirectCall
#       likelihood accumulation regression through a MultiSampler.
# ======================================================================

# Helper chooser that routes every clock to a single sub-sampler (key 1). Used
# by the MultiSampler routing regressions below.
module CNRMultiHelp
using CompetingClocks
using CompetingClocks: SamplerChoice, choose_sampler
using Distributions: UnivariateDistribution
struct AllToOne <: SamplerChoice{Int64,Int64} end
function CompetingClocks.choose_sampler(
    ::AllToOne, clock::Int64, distribution::UnivariateDistribution
)::Int64
    return 1
end
end


@safetestset CNR_next_purity = "CNR next() is pure and consumes no RNG" begin
    using CompetingClocks: CombinedNextReaction, enable!, next
    using Random: Xoshiro
    using Distributions: Weibull, Exponential

    rng = Xoshiro(11)
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)
    enable!(sampler, 2, Exponential(1.0), 0.0, 0.0)

    # Snapshot the RNG, snapshot the sampler's internal state, then call next
    # several times with no intervening enable!/disable!/fire!.
    entry1_before = sampler.transition_entry[1]
    entry2_before = sampler.transition_entry[2]
    # The sampler owns its randomness now; a draw shows up as a bumped per-key
    # occurrence count in its keyed streams, so we observe consumption there
    # rather than by watching an externally-passed rng advance.
    counts_before = copy(sampler.streams.counts)

    r1 = next(sampler, 0.0)
    r2 = next(sampler, 0.0)
    r3 = next(sampler, 0.0)

    # Purity: repeated next() returns an identical (when, key) reservation.
    @test r1 == r2
    @test r2 == r3

    # next() mutated no stored survival (the old code zeroed the returned clock).
    @test sampler.transition_entry[1] == entry1_before
    @test sampler.transition_entry[2] == entry2_before

    # next() consumed no random numbers: no per-key stream advanced.
    @test sampler.streams.counts == counts_before
end


@safetestset CNR_fired_path_redraws = "CNR fire! consumes draw; re-enable redraws fresh" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, fire!
    using CompetingClocks: get_survival_zero, sampling_space
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(22)
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0)

    when, which = next(sampler, 0.0)
    @test which == 1

    fire!(sampler, 1, when)
    # fire! removed the clock from the queue and zeroed its survival completely.
    rec = sampler.transition_entry[1]
    @test rec.heap_handle == 0
    @test rec.survival == get_survival_zero(sampling_space(Exponential))

    # Re-enabling (with any distribution) must take the fresh-redraw branch,
    # which draws from clock 1's own keyed stream (bumping its occurrence count).
    c0 = get(sampler.streams.counts, 1, 0)
    enable!(sampler, 1, Exponential(2.0), when, when)
    @test get(sampler.streams.counts, 1, 0) == c0 + 1   # stream advanced => fresh redraw
end


@safetestset CNR_loser_path_reuses = "CNR loser (returned-but-not-fired) re-enable reuses draw" begin
    # This is the bug the change fixes. Before the fix, next() marked the
    # returned clock's survival to zero, so re-enabling that (non-fired) clock
    # took the fresh-redraw branch, destroying draw reuse. After the fix, next()
    # is pure and the re-enable REUSES the draw via consume_survival +
    # sample_by_inversion, consuming no RNG.
    using CompetingClocks: CombinedNextReaction, enable!, next
    using CompetingClocks: consume_survival, sample_by_inversion, sampling_space, LogSampling
    using Random: Xoshiro
    using Distributions: Weibull

    rng = Xoshiro(33)
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)
    enable!(sampler, 2, Weibull(1.5, 2.0), 0.0, 0.0)

    when, loser = next(sampler, 0.0)   # the returned (minimum) clock; do NOT fire it

    # Re-enable that same, non-fired clock with a CHANGED distribution at a later
    # time, keeping te the same so we hit the consume_survival reuse branch.
    newdist = Weibull(3.0, 1.5)
    S = sampling_space(newdist)             # LogSampling
    @test S === LogSampling
    when_re = 0.5
    rec = sampler.transition_entry[loser]
    te = rec.te

    # Compute the expected reused firing time analytically with the package's own
    # cache functions.
    survival_remain = consume_survival(rec, rec.distribution, S, when_re)
    expected_tau = sample_by_inversion(newdist, S, te, when_re, survival_remain)

    c0 = get(sampler.streams.counts, loser, 0)
    enable!(sampler, loser, newdist, te, when_re)

    # Reuse => no draw from the clock's stream, and the stored firing time matches
    # the analytic inversion of the recorded survival.
    @test get(sampler.streams.counts, loser, 0) == c0
    @test sampler[loser] ≈ expected_tau
end


@safetestset CNR_disable_without_fire_reuses = "CNR disable (no fire) preserves survival for reuse" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, disable!
    using CompetingClocks: sample_by_inversion, sampling_space, get_survival_zero, LogSampling
    using Random: Xoshiro
    using Distributions: Weibull

    rng = Xoshiro(44)
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)
    next(sampler, 0.0)          # queried but not fired
    disable!(sampler, 1, 0.3)        # preserves remaining survival, heap_handle -> 0

    rec = sampler.transition_entry[1]
    @test rec.heap_handle == 0
    @test rec.survival > get_survival_zero(LogSampling)   # remaining survival kept (> -Inf)

    newdist = Weibull(2.0, 1.0)
    S = sampling_space(newdist)
    when_re = 0.5
    # Disabled-branch reuse uses record.survival directly (no consume_survival).
    expected_tau = sample_by_inversion(newdist, S, rec.te, when_re, rec.survival)

    c0 = get(sampler.streams.counts, 1, 0)
    enable!(sampler, 1, newdist, rec.te, when_re)

    @test get(sampler.streams.counts, 1, 0) == c0   # reuse => no draw
    @test sampler[1] ≈ expected_tau
end


@safetestset CNR_consume_survival_analytic = "CNR consume_survival analytic in both spaces" begin
    using CompetingClocks: NRTransition, consume_survival, LinearSampling, LogSampling
    using Distributions: Uniform, Weibull, ccdf, logccdf

    te = 0.0
    t0 = 1.0
    tn = 3.0

    # --- LinearSampling (Gibson-Bruck): survival * (S(te,t0) / S(te,tn)) ---
    dist_lin = Uniform(0.0, 10.0)
    s0 = 0.4
    rec_lin = NRTransition{Float64}(1, s0, dist_lin, te, t0)
    expected_lin = s0 * (ccdf(dist_lin, t0 - te) / ccdf(dist_lin, tn - te))
    @test consume_survival(rec_lin, dist_lin, LinearSampling, tn) ≈ expected_lin

    # --- LogSampling (Anderson): survival - (logS(te,tn) - logS(te,t0)) ---
    dist_log = Weibull(2.0, 1.0)
    sl0 = -0.7
    rec_log = NRTransition{Float64}(1, sl0, dist_log, te, t0)
    expected_log = sl0 - (logccdf(dist_log, tn - te) - logccdf(dist_log, t0 - te))
    @test consume_survival(rec_log, dist_log, LogSampling, tn) ≈ expected_log

    # --- Edge case: te at/after t0 and tn -> the conditional survivals are 1
    #     (linear) / 0 (log), so consume_survival returns the stored survival. ---
    te2 = 5.0
    rec_lin2 = NRTransition{Float64}(1, s0, dist_lin, te2, t0)
    @test consume_survival(rec_lin2, dist_lin, LinearSampling, tn) ≈ s0
    rec_log2 = NRTransition{Float64}(1, sl0, dist_log, te2, t0)
    @test consume_survival(rec_log2, dist_log, LogSampling, tn) ≈ sl0
end


@safetestset CNR_same_te_same_dist_noop = "CNR re-enable with same te and dist is a no-op" begin
    using CompetingClocks: CombinedNextReaction, enable!, next
    using Random: Xoshiro
    using Distributions: Weibull

    rng = Xoshiro(55)
    sampler = CombinedNextReaction{Int,Float64}()
    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)

    t_before = sampler[1]
    rec_before = sampler.transition_entry[1]
    c0 = get(sampler.streams.counts, 1, 0)

    # Re-enabling an already-enabled clock with the identical distribution and
    # enabling time must change nothing and consume no randomness.
    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)

    @test sampler[1] == t_before
    rec_after = sampler.transition_entry[1]
    @test rec_after.survival == rec_before.survival
    @test rec_after.heap_handle == rec_before.heap_handle
    @test get(sampler.streams.counts, 1, 0) == c0
end


@safetestset CNR_multisampler_fire_routing = "MultiSampler fire! reaches CNR sub-sampler" begin
    # Regression: MultiSampler.fire! must forward to the sub-sampler's fire!
    # (not disable!). Proof: after next() (which no longer mutates), fire!
    # through the MultiSampler must fully consume the CNR draw, so re-enabling
    # the fired clock redraws fresh (consumes RNG). Had it routed to disable!,
    # the remaining survival would be reused and no RNG consumed.
    using CompetingClocks: CombinedNextReaction, MultiSampler, enable!, next, fire!
    using CompetingClocks: get_survival_zero, LogSampling
    using ..CNRMultiHelp: AllToOne
    using Random: Xoshiro
    using Distributions: Weibull

    rng = Xoshiro(66)
    sampler = MultiSampler{Int64,Int64,Float64}(AllToOne())
    sampler[1] = CombinedNextReaction{Int64,Float64}()

    enable!(sampler, 1, Weibull(2.0, 1.0), 0.0, 0.0)
    when, which = next(sampler, 0.0)
    @test which == 1

    fire!(sampler, 1, when)

    sub = sampler.propagator[1]
    @test sub.transition_entry[1].survival == get_survival_zero(LogSampling)

    c0 = get(sub.streams.counts, 1, 0)
    enable!(sampler, 1, Weibull(2.0, 1.0), when, when)
    # fresh redraw (a bumped stream count on the sub-sampler) => fire! reached it
    @test get(sub.streams.counts, 1, 0) == c0 + 1
end


@safetestset CNR_multisampler_directcall_likelihood = "MultiSampler fire! accumulates DirectCall likelihood" begin
    # Regression for the separate DirectCall bug: DirectCall.fire! accumulates
    # sampler-native step log-likelihood beyond disable!. MultiSampler.fire!
    # previously routed to disable! and silently skipped it. Now fire! reaches
    # the sub-sampler, so the accumulated likelihood matches firing a bare
    # DirectCall directly.
    using CompetingClocks: DirectCall, MultiSampler, enable!, fire!
    using ..CNRMultiHelp: AllToOne
    using Random: Xoshiro
    using Distributions: Exponential

    # Bare DirectCall reference.
    rng1 = Xoshiro(77)
    dc = DirectCall{Int64,Float64}(trajectory=true)
    enable!(dc, 1, Exponential(1.0), 0.0, 0.0)
    enable!(dc, 2, Exponential(2.0), 0.0, 0.0)
    fire!(dc, 1, 0.5)
    bare_ll = dc.log_likelihood

    # Same firing, but routed through a MultiSampler.
    rng2 = Xoshiro(77)
    ms = MultiSampler{Int64,Int64,Float64}(AllToOne())
    ms[1] = DirectCall{Int64,Float64}(trajectory=true)
    enable!(ms, 1, Exponential(1.0), 0.0, 0.0)
    enable!(ms, 2, Exponential(2.0), 0.0, 0.0)
    fire!(ms, 1, 0.5)
    ms_ll = ms.propagator[1].log_likelihood

    @test ms_ll == bare_ll
    @test ms_ll != 0.0    # likelihood really was accumulated (disable! would leave it 0)
end


@safetestset CNR_steploglikelihood_enabled_only = "CNR steploglikelihood excludes fired and disabled clocks" begin
    # Fired and disabled clocks are retained in transition_entry (heap_handle == 0)
    # for draw reuse, but they must not contribute spurious survival terms to the
    # step log-likelihood. The TrackWatcher, which only tracks enabled clocks,
    # provides the ground-truth value.
    using CompetingClocks: CombinedNextReaction, TrackWatcher
    using CompetingClocks: enable!, fire!, disable!, steploglikelihood
    using Distributions: Exponential, logpdf, logccdf
    using Random: Xoshiro

    # Analytic value for the exponential case: the survivor of the disabled/fired
    # clock 1 must not appear, leaving only clock 2's contribution.
    analytic = logpdf(Exponential(1.0), 1.0) - logccdf(Exponential(1.0), 0.5)
    @test analytic ≈ -0.5

    # --- Case 1: fire! retains clock 1 with zeroed survival. ---
    rng = Xoshiro(101)
    nr = CombinedNextReaction{Int,Float64}()
    tw = TrackWatcher{Int,Float64}()
    for dst in (nr, tw)
        enable!(dst, 1, Exponential(1.0), 0.0, 0.0)
        enable!(dst, 2, Exponential(1.0), 0.0, 0.0)
        fire!(dst, 1, 0.5)
    end
    @test steploglikelihood(nr, 0.5, 1.0, 2) ≈ steploglikelihood(tw, 0.5, 1.0, 2)
    @test steploglikelihood(nr, 0.5, 1.0, 2) ≈ analytic
    @test steploglikelihood(tw, 0.5, 1.0, 2) ≈ analytic

    # --- Case 2: disable! retains clock 1 with nonzero remaining survival. ---
    rng = Xoshiro(202)
    nr = CombinedNextReaction{Int,Float64}()
    tw = TrackWatcher{Int,Float64}()
    for dst in (nr, tw)
        enable!(dst, 1, Exponential(1.0), 0.0, 0.0)
        enable!(dst, 2, Exponential(1.0), 0.0, 0.0)
        disable!(dst, 1, 0.5)
    end
    @test steploglikelihood(nr, 0.5, 1.0, 2) ≈ steploglikelihood(tw, 0.5, 1.0, 2)
    @test steploglikelihood(nr, 0.5, 1.0, 2) ≈ analytic

    # --- Case 3: re-enabling clock 1 makes it participate again. ---
    enable!(nr, 1, Exponential(1.0), 0.5, 0.5)
    enable!(tw, 1, Exponential(1.0), 0.5, 0.5)
    @test steploglikelihood(nr, 0.5, 1.0, 2) ≈ steploglikelihood(tw, 0.5, 1.0, 2)
    # Now clock 1 contributes survival again, so the value drops below analytic.
    @test steploglikelihood(nr, 0.5, 1.0, 2) < analytic
end


@safetestset combinednr_jitter_with_disabled_entry =
    "combinednr: jitter! on a sampler holding a retained-disabled entry does not throw and extinguishes the banked survival so a re-enable draws fresh" begin
    using CompetingClocks
    using CompetingClocks: CombinedNextReaction, jitter!
    using Distributions: Weibull

    # A disabled (not fired) clock banks residual survival for reuse
    # (heap_handle == 0 entry kept in the table). jitter!'s purpose is
    # decorrelation, which must extinguish that bank too; NRTransition is
    # immutable, so the old in-place assignment threw a setfield! error the
    # moment any disabled entry existed.
    nr = CombinedNextReaction{Int,Float64}(424242)
    enable!(nr, 1, Weibull(1.7, 2.0), 0.0, 0.0)
    enable!(nr, 2, Weibull(1.3, 1.5), 0.0, 0.0)
    disable!(nr, 1, 0.4)                    # banks clock 1's residual survival
    jitter!(nr, 0.4)                        # threw before the fix
    # The bank is gone: re-enabling clock 1 draws a fresh schedule rather than
    # reusing pre-jitter randomness, and the sampler stays fully operational.
    enable!(nr, 1, Weibull(1.7, 2.0), 0.4, 0.4)
    when, which = next(nr, 0.4)
    @test which in (1, 2)
    @test when > 0.4
end
