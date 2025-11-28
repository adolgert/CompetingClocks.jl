using SafeTestsets


@safetestset CombinedNextReactionSmoke = "CombinedNextReaction reaction does basic things" begin
    using Distributions
    using Random
    using CompetingClocks: CombinedNextReaction, next, enable!, disable!, reset!

    rng = MersenneTwister(349827)
    for i in 1:100
        sampler = CombinedNextReaction{String,Float64}()
        @test next(sampler, 3.0, rng)[2] === nothing
        enable!(sampler, "walk home", Exponential(1.5), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] == "walk home"
        enable!(sampler, "run", Gamma(1, 3), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        enable!(sampler, "walk to sandwich shop", Weibull(2, 1), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run", "walk to sandwich shop"]
        disable!(sampler, "walk to sandwich shop", 1.7)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        reset!(sampler)
    end
end

@safetestset CombinedNextReaction_interface = "CombinedNextReaction basic interface" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    sampler = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[1]
    @test keytype(sampler) <: Int64

    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[1] == 7.9

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "1")

    disable!(sampler, 1, 0.0)

    @test_throws BoundsError sampler[1]
    @test sampler[2] == 12.3

end


@safetestset CombinedNextReaction_pair_keys = "CombinedNextReaction pair keys" begin
    using CompetingClocks
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
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[:S=>1] == 7.9

    @test haskey(sampler, :S => 1)
    @test !haskey(sampler, :I => 17)

    disable!(sampler, :S => 1, 0.0)

    @test_throws BoundsError sampler[:S=>1]
    @test sampler[:S=>2] == 12.3

end

@safetestset CombinedNextReaction_copy = "CombinedNextReaction copy" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    dst = clone(src)
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0, rng)
    enable!(src, 38, Exponential(), 0.0, 0.0, rng)
    enable!(dst, 29, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(src) == 2
    @test length(dst) == 2
    # Changing src doesn't change dst.
    enable!(src, 48, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    # Changing dst doesn't change src.
    enable!(dst, 49, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3
end


@safetestset CombinedNextReaction_set = "CombinedNextReaction set" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0, rng)
    enable!(src, 38, Exponential(), 0.0, 0.0, rng)
    enable!(src, 29, Exponential(), 0.0, 0.0, rng)
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
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)

    # Get current next event
    t1, k1 = next(sampler, 0.0, rng)

    # Jitter should resample all clocks
    jitter!(sampler, 0.5, rng)

    # Get new next event - times should be different after jitter
    t2, k2 = next(sampler, 0.5, rng)
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
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.5, rng)
    t, k = next(sampler, 0.5, rng)
    @test t >= 0.5  # Fire time must be >= current time (truncated dist guarantees this)

    # Test with LogSampling distribution (Gamma)
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.5, rng)
    t2, k2 = next(sampler2, 0.5, rng)
    @test t2 >= 0.5

    # Test with Weibull (also LogSampling)
    sampler3 = CombinedNextReaction{Int,Float64}()
    enable!(sampler3, 1, Weibull(2.0, 1.0), 0.0, 0.5, rng)
    t3, k3 = next(sampler3, 0.5, rng)
    @test t3 >= 0.5
end


@safetestset CombinedNextReaction_reenable = "CombinedNextReaction re-enable after fire" begin
    using CompetingClocks: CombinedNextReaction, enable!, next, fire!, disable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable and fire a clock
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    t, k = next(sampler, 0.0, rng)
    fire!(sampler, 1, t)

    # Re-enable the same clock (should exercise heap_handle > 0 but survival == 0 branch)
    enable!(sampler, 1, Exponential(1.0), t, t, rng)
    t2, k2 = next(sampler, t, rng)
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
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    t1, k1 = next(sampler, 0.0, rng)

    # Update the same clock with different distribution (exercises consume_survival)
    enable!(sampler, 1, Exponential(2.0), 0.0, 0.5, rng)
    t2, k2 = next(sampler, 0.5, rng)
    @test t2 >= 0.5

    # Update with a LogSampling distribution to exercise that branch
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.0, rng)
    t3, k3 = next(sampler2, 0.0, rng)
    enable!(sampler2, 1, Gamma(3.0, 1.0), 0.0, 0.5, rng)
    t4, k4 = next(sampler2, 0.5, rng)
    @test t4 >= 0.5
end


@safetestset CombinedNextReaction_isenabled = "CombinedNextReaction isenabled" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, isenabled
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(678901)
    sampler = CombinedNextReaction{Int,Float64}()

    @test !isenabled(sampler, 1)

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
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
    enable!(sampler1, 1, truncated(Normal(5.0, 1.0), 0, Inf), 0.0, 0.0, rng)
    t1, k1 = next(sampler1, 0.0, rng)
    @test t1 > 0

    # Test with Uniform distribution
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Uniform(1.0, 5.0), 0.0, 0.0, rng)
    t2, k2 = next(sampler2, 0.0, rng)
    @test 1.0 <= t2 <= 5.0

    # Test with LogNormal
    sampler3 = CombinedNextReaction{Int,Float64}()
    enable!(sampler3, 1, LogNormal(0.0, 1.0), 0.0, 0.0, rng)
    t3, k3 = next(sampler3, 0.0, rng)
    @test t3 > 0

    # Test with Pareto
    sampler4 = CombinedNextReaction{Int,Float64}()
    enable!(sampler4, 1, Pareto(1.0), 0.0, 0.0, rng)
    t4, k4 = next(sampler4, 0.0, rng)
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
    enable!(sampler, 1, Exponential(1.0), 2.0, 0.0, rng)
    t, k = next(sampler, 0.0, rng)
    @test t >= 2.0  # Fire time must be after the enabling time
end


@safetestset CombinedNextReaction_disable_reenable = "CombinedNextReaction disable then re-enable" begin
    using CompetingClocks: CombinedNextReaction, enable!, disable!, next
    using Random: Xoshiro
    using Distributions: Exponential, Gamma

    rng = Xoshiro(901234)
    sampler = CombinedNextReaction{Int,Float64}()

    # Enable, disable, then re-enable (should use remaining survival)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    disable!(sampler, 1, 0.3)

    # Re-enable the disabled clock (heap_handle == 0 but has remaining survival)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.5, rng)
    t, k = next(sampler, 0.5, rng)
    @test k == 1
    @test t >= 0.5

    # Same test with LogSampling distribution
    sampler2 = CombinedNextReaction{Int,Float64}()
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.0, rng)
    disable!(sampler2, 1, 0.3)
    enable!(sampler2, 1, Gamma(2.0, 1.0), 0.0, 0.5, rng)
    t2, k2 = next(sampler2, 0.5, rng)
    @test k2 == 1
    @test t2 >= 0.5
end