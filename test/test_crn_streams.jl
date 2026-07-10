using SafeTestsets

# The successor to the old record/replay CommonRandom machinery. Because the
# sampler owns per-clock keyed streams, common random numbers need no recorder
# and no replay mode: two samplers built from the SAME seed draw identically for
# the same clock, whatever order their events fire in. These tests are the
# evidence for that design claim — an adversarial interleaving that still couples
# per clock, and a coupled θ / θ+h pair that reduces variance against independent
# seeds on a machine-repair mean.


@safetestset crn_streams_order_independence =
    "two samplers from one seed, driven through different event orders, fire each clock at the identical time where its enabled history matches" begin
    using CompetingClocks: FirstToFire, enable!, disable!, next, fire!
    using Distributions: Weibull, Exponential, Gamma

    dists = Dict(:a => Weibull(1.5, 1.0), :b => Exponential(1.3),
                 :c => Gamma(2.0, 0.7), :d => Weibull(2.0, 0.9))

    # Run 1 enables in one order; Run 2 in the reverse order. Each clock is drawn
    # once, at enable, from ITS OWN stream, so its scheduled firing time cannot
    # depend on when the other clocks were enabled.
    s1 = FirstToFire{Symbol,Float64}(0xC0FFEE)
    for k in (:a, :b, :c, :d)
        enable!(s1, k, dists[k], 0.0, 0.0)
    end

    s2 = FirstToFire{Symbol,Float64}(0xC0FFEE)
    for k in (:d, :c, :b, :a)   # reverse enable order
        enable!(s2, k, dists[k], 0.0, 0.0)
    end

    for k in (:a, :b, :c, :d)
        @test s1[k] == s2[k]
    end

    # And the full race sequences agree, because identical per-clock times imply
    # an identical ordering.
    seq1 = Tuple{Float64,Symbol}[]; seq2 = Tuple{Float64,Symbol}[]
    for _ in 1:4
        w1, k1 = next(s1, 0.0); push!(seq1, (w1, k1)); fire!(s1, k1, w1)
        w2, k2 = next(s2, 0.0); push!(seq2, (w2, k2)); fire!(s2, k2, w2)
    end
    @test seq1 == seq2
end


@safetestset crn_streams_reuse_across_interleave =
    "a clock disabled and re-enabled in different global orders keeps a coupled per-clock draw on CombinedNextReaction" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, next
    using Distributions: Weibull

    # Same seed, same per-clock history for :x (enable, disable at 0.3, re-enable
    # at 0.5), but the runs touch :y and :z in different orders in between.
    s1 = CombinedNextReaction{Symbol,Float64}(0x5EED5)
    enable!(s1, :x, Weibull(1.7, 1.0), 0.0, 0.0)
    enable!(s1, :y, Weibull(1.2, 2.0), 0.0, 0.0)
    disable!(s1, :x, 0.3)
    enable!(s1, :z, Weibull(2.0, 1.0), 0.0, 0.0)
    enable!(s1, :x, Weibull(1.7, 1.0), 0.0, 0.5)

    s2 = CombinedNextReaction{Symbol,Float64}(0x5EED5)
    enable!(s2, :z, Weibull(2.0, 1.0), 0.0, 0.0)   # different order for the others
    enable!(s2, :y, Weibull(1.2, 2.0), 0.0, 0.0)
    enable!(s2, :x, Weibull(1.7, 1.0), 0.0, 0.0)
    disable!(s2, :x, 0.3)
    enable!(s2, :x, Weibull(1.7, 1.0), 0.0, 0.5)

    @test s1[:x] == s2[:x]
end


# ---------------------------------------------------------------------------
# Coupled θ / θ+h experiment on a machine-repair model. With CRN (same seed for
# the θ and θ+h runs) each machine's k-th break/repair draw is the SAME underlying
# variate in both runs, so the paired difference has far less variance than under
# independent seeds. This is the variance reduction the old CommonRandom aimed at,
# now falling out of stream ownership.
# ---------------------------------------------------------------------------
@safetestset crn_streams_variance_reduction =
    "coupling a machine-repair mean at theta and theta+h from one seed reduces the paired-difference variance far below independent seeds" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, next, fire!
    using Distributions: Exponential
    using Statistics: mean, var, std

    N = 5
    T = 50.0
    MU = 0.5           # repair-time scale
    THETA = 1.0        # break-time scale
    H = 0.05

    # Count repair completions by horizon T. break-time scale is theta; a shorter
    # scale => more break/repair cycles => more repairs, so f is theta-sensitive.
    function repairs_by_horizon(seed, theta)
        s = CombinedNextReaction{Tuple{Int,Symbol},Float64}(seed)
        for i in 1:N
            enable!(s, (i, :break), Exponential(theta), 0.0, 0.0)
        end
        t = 0.0
        repairs = 0
        while true
            when, which = next(s, t)
            (which === nothing || when > T) && break
            i, kind = which
            fire!(s, which, when)
            if kind === :break
                enable!(s, (i, :repair), Exponential(MU), when, when)
            else
                repairs += 1
                enable!(s, (i, :break), Exponential(theta), when, when)
            end
            t = when
        end
        return repairs
    end

    reps = 3000
    crn = Vector{Float64}(undef, reps)
    indep = Vector{Float64}(undef, reps)
    for r in 1:reps
        # CRN: the SAME seed drives both theta and theta+h.
        crn[r] = repairs_by_horizon(r, THETA + H) - repairs_by_horizon(r, THETA)
        # Independent: disjoint seeds for the two runs.
        indep[r] = repairs_by_horizon(10_000_000 + r, THETA + H) -
                   repairs_by_horizon(20_000_000 + r, THETA)
    end

    var_crn = var(crn)
    var_indep = var(indep)
    @test var_crn > 0                       # coupling is not degenerate/exact
    @test var_indep / var_crn > 5.0         # substantial variance reduction

    # Both estimate the same mean difference; they must agree within a few SE.
    se = hypot(std(crn) / sqrt(reps), std(indep) / sqrt(reps))
    @test abs(mean(crn) - mean(indep)) < 4 * se
    # And the mean effect is resolved well enough by the CRN estimator to matter.
    @test abs(mean(crn)) > 3 * (std(crn) / sqrt(reps))
end
