using SafeTestsets


@safetestset firsttofire_smoke = "FirstToFire smoke" begin
    using CompetingClocks: FirstToFire, enable!, disable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    sampler = FirstToFire{Int64,Float64}()
    rng = Xoshiro(90422342)
    enabled = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0, rng)
        push!(enabled, clock_id)
    end
    when, which = next(sampler, 0.0, rng)
    disable!(sampler, which, when)
    delete!(enabled, which)
    when, which = next(sampler, when, rng)
    @test when > 0.0
    @test 1 <= which
    @test which <= 5
    @test which ∈ enabled
end


@safetestset FirstToFire_interface = "FirstToFire basic interface" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    rng = Xoshiro(123)

    propagator = FirstToFire{Int64,Float64}()

    @test length(propagator) == 0
    @test length(keys(propagator)) == 0
    @test_throws KeyError propagator[1]
    @test keytype(propagator) <: Int64

    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        enable!(propagator, clock, Dirac(when_fire), 0.0, 0.0, rng)
    end

    @test length(propagator) == 5
    @test length(keys(propagator)) == 5
    @test propagator[1] == 7.9

    @test haskey(propagator, 1)
    @test !haskey(propagator, 1_000)
    @test !haskey(propagator, "1")

    disable!(propagator, 1, 0.0)

    @test_throws KeyError propagator[1]
    @test propagator[2] == 12.3

end


@safetestset FirstToFire_copy = "FirstToFire copy" begin
    using CompetingClocks: FirstToFire, enable!, next, copy_clocks!
    using Random: MersenneTwister
    using Distributions: Exponential

    src = FirstToFire{Int,Float64}()
    dst = FirstToFire{Int,Float64}()
    rng = MersenneTwister(90422342)
    enable!(src, 1, Exponential(), 0.0, 0.0, rng)
    enable!(src, 2, Exponential(), 0.0, 0.0, rng)
    enable!(dst, 3, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(dst) == 2
    enable!(src, 5, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    enable!(dst, 6, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3
end


@safetestset FirstToFire_clone = "FirstToFire clone" begin
    using CompetingClocks: FirstToFire, enable!, clone
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(234567)
    sampler = FirstToFire{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)

    cloned = clone(sampler)
    @test length(cloned) == 0  # cloned is empty
    @test length(sampler) == 2  # original unchanged
end


@safetestset FirstToFire_jitter = "FirstToFire jitter!" begin
    using CompetingClocks: FirstToFire, enable!, next, jitter!
    using Random: Xoshiro
    using Distributions: Exponential, Gamma

    rng = Xoshiro(345678)
    sampler = FirstToFire{Int,Float64}()

    # Enable clocks
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Gamma(2.0, 1.0), 0.0, 0.0, rng)

    # Get initial next event
    t1, k1 = next(sampler, 0.0, rng)

    # Jitter resamples all clocks - times should change
    jitter!(sampler, 0.5, rng)

    # Get new next event after jitter
    t2, k2 = next(sampler, 0.5, rng)
    @test t2 >= 0.5  # New times are after jitter time

    # Test jitter with te < when (truncated distribution branch)
    sampler2 = FirstToFire{Int,Float64}()
    enable!(sampler2, 1, Exponential(1.0), 0.0, 0.0, rng)
    jitter!(sampler2, 0.3, rng)  # when > te, so uses truncated distribution
    t3, k3 = next(sampler2, 0.3, rng)
    @test t3 >= 0.3

    # Test jitter with te >= when (non-truncated distribution branch)
    sampler3 = FirstToFire{Int,Float64}()
    enable!(sampler3, 1, Exponential(1.0), 1.0, 1.0, rng)  # te = 1.0
    jitter!(sampler3, 0.5, rng)  # when < te, so uses regular distribution
    t4, k4 = next(sampler3, 0.5, rng)
    @test t4 >= 1.0  # Fire time should be >= te
end


@safetestset FirstToFire_haskey_wrong_type = "FirstToFire haskey with wrong type" begin
    using CompetingClocks: FirstToFire, enable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = FirstToFire{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)

    # haskey with correct type
    @test haskey(sampler, 1) == true
    @test haskey(sampler, 999) == false

    # haskey with wrong type should return false (not throw)
    @test haskey(sampler, "wrong_type") == false
    @test haskey(sampler, :symbol) == false
    @test haskey(sampler, 1.5) == false
end
