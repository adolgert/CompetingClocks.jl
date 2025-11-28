using SafeTestsets


@safetestset rssa_smoke = "RSSA smoke" begin
    using CompetingClocks: RSSA, enable!, disable!, next
    using Random: MersenneTwister
    using Distributions: Exponential

    rng = MersenneTwister(90422342)
    min_when = Inf
    seen = Set{Int}()
    sample_time = 0.5
    for i in 1:100
        sampler = RSSA{Int,Float64}()
        enable!(sampler, 1, Exponential(1.7), 0.0, 0.0, rng)
        enable!(sampler, 2, Exponential(4.5), 0.0, 0.0, rng)
        enable!(sampler, 3, Exponential(2.0), 0.0, 0.0, rng)
        disable!(sampler, 2, sample_time)
        when, which = next(sampler, sample_time, rng)
        push!(seen, which)
        @test when > sample_time
        min_when = min(min_when, when)
    end
    @test min_when < sample_time + 0.03
    @test seen == Set([1, 3])
end


@safetestset rssa_interface = "RSSA basic interface" begin
    using CompetingClocks: RSSA, enable!, disable!, next, reset!, enabled, isenabled
    using Random: Xoshiro
    using Distributions: Exponential

    sampler = RSSA{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test keytype(sampler) <: Int64

    for (clock, rate) in [(1, 0.5), (2, 1.0), (3, 2.0), (4, 10.0), (5, 0.1)]
        enable!(sampler, clock, Exponential(1/rate), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test length(enabled(sampler)) == 5

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "wrong_type")
    @test isenabled(sampler, 1)
    @test !isenabled(sampler, 1_000)

    disable!(sampler, 1, 0.0)
    @test !haskey(sampler, 1)
    @test !isenabled(sampler, 1)
    reset!(sampler)
    @test length(sampler) == 0
end


@safetestset rssa_empty = "RSSA empty" begin
    using CompetingClocks: RSSA, next
    using Random: MersenneTwister

    rng = MersenneTwister(90422342)
    sampler = RSSA{Int,Float64}()
    when, which = next(sampler, 5.7, rng)
    @test when == Inf
    @test which === nothing
end


@safetestset rssa_clone = "RSSA clone" begin
    using CompetingClocks: RSSA, enable!, clone
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(234567)
    sampler = RSSA{Int,Float64}(bound_factor=1.1)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)
    cloned = clone(sampler)
    @test cloned.bound_factor == 1.1
    @test length(cloned) == 0  # cloned is empty
    @test length(sampler) == 2  # original unchanged
end


@safetestset rssa_copy = "RSSA copy_clocks!" begin
    using CompetingClocks: RSSA, enable!, next, copy_clocks!
    using Random: MersenneTwister
    using Distributions: Exponential

    src = RSSA{Int,Float64}(bound_factor=1.2)
    dst = RSSA{Int,Float64}(bound_factor=1.3)
    rng = MersenneTwister(90422342)
    enable!(src, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(src, 2, Exponential(1.0), 0.0, 0.0, rng)
    enable!(dst, 3, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(dst) == 2
    @test dst.bound_factor == 1.2  # copied from src
    enable!(src, 5, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    enable!(dst, 6, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3
end


@safetestset rssa_set_bound = "RSSA set_bound!" begin
    using CompetingClocks: RSSA, enable!, set_bound!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(345678)
    sampler = RSSA{Int,Float64}(bound_factor=1.0)  # No inflation

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0
    @test sampler.a[1] == 1.0
    @test sampler.abar[1] == 1.0  # bound_factor=1.0, so abar = rate

    # Increase the bound
    set_bound!(sampler, 1, 2.0)
    @test sampler.abar[1] == 2.0
    @test sampler.Abar == 2.0

    # Try to set bound below true rate - should clamp to rate
    set_bound!(sampler, 1, 0.5)
    @test sampler.abar[1] == 1.0  # clamped to a[1]

    # Test KeyError for non-existent clock
    @test_throws KeyError set_bound!(sampler, 999, 1.0)
end


@safetestset rssa_set_global_bound_factor = "RSSA set_global_bound_factor!" begin
    using CompetingClocks: RSSA, enable!, set_global_bound_factor!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = RSSA{Int,Float64}(bound_factor=1.0)

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0
    enable!(sampler, 2, Exponential(0.5), 0.0, 0.0, rng)  # rate = 2.0
    @test sampler.abar[1] == 1.0
    @test sampler.abar[2] == 2.0

    # Change global bound factor
    set_global_bound_factor!(sampler, 1.5)
    @test sampler.bound_factor == 1.5
    @test sampler.abar[1] == 1.5  # 1.0 * 1.5
    @test sampler.abar[2] == 3.0  # 2.0 * 1.5
    @test sampler.Abar == 4.5

    # Try to set bound_factor < 1 - should clamp to 1.0
    set_global_bound_factor!(sampler, 0.5)
    @test sampler.bound_factor == 1.0
end


@safetestset rssa_jitter = "RSSA jitter!" begin
    using CompetingClocks: RSSA, enable!, next, jitter!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(567890)
    sampler = RSSA{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)

    # Get next event to populate cache
    t1, k1 = next(sampler, 0.0, rng)
    @test sampler.cached_next !== nothing

    # Jitter invalidates cache
    jitter!(sampler, 0.0, rng)
    @test sampler.cached_next === nothing

    # Next call should give a new sample
    t2, k2 = next(sampler, 0.0, rng)
    @test sampler.cached_next !== nothing
end


@safetestset rssa_non_exponential = "RSSA non-exponential error" begin
    using CompetingClocks: RSSA, enable!
    using Random: Xoshiro
    using Distributions: Gamma, Weibull

    rng = Xoshiro(678901)
    sampler = RSSA{Int,Float64}()

    @test_throws ArgumentError enable!(sampler, 1, Gamma(1.0, 1.0), 0.0, 0.0, rng)
    @test_throws ArgumentError enable!(sampler, 2, Weibull(1.0, 1.0), 0.0, 0.0, rng)
end


@safetestset rssa_rate_update = "RSSA rate update" begin
    using CompetingClocks: RSSA, enable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(789012)
    sampler = RSSA{Int,Float64}(bound_factor=1.0)

    # Enable with rate 1.0
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0
    @test sampler.a[1] == 1.0
    @test sampler.abar[1] == 1.0

    # Update to higher rate - bound should increase
    enable!(sampler, 1, Exponential(0.5), 0.0, 0.0, rng)  # rate = 2.0
    @test sampler.a[1] == 2.0
    @test sampler.abar[1] == 2.0  # bound increases to match rate

    # Update to lower rate - bound stays at previous level
    enable!(sampler, 1, Exponential(2.0), 0.0, 0.0, rng)  # rate = 0.5
    @test sampler.a[1] == 0.5
    @test sampler.abar[1] == 2.0  # bound doesn't decrease automatically
end


@safetestset rssa_getindex = "RSSA getindex" begin
    using CompetingClocks: RSSA, enable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(890123)
    sampler = RSSA{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)

    # No cached_next yet - should throw KeyError
    @test_throws KeyError sampler[1]

    # After next(), cached_next is set
    t, k = next(sampler, 0.0, rng)
    @test k == 1  # only one clock enabled
    @test sampler[1] == t

    # Add another clock, ask for non-cached clock
    enable!(sampler, 2, Exponential(1.0), 0.0, 0.0, rng)
    t2, k2 = next(sampler, 0.0, rng)
    # Only the cached clock should be accessible
    @test sampler[k2] == t2
    if k2 == 1
        @test_throws KeyError sampler[2]
    else
        @test_throws KeyError sampler[1]
    end
end


@safetestset rssa_cached_next = "RSSA cached next idempotent" begin
    using CompetingClocks: RSSA, enable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(901234)
    sampler = RSSA{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)

    # First call sets the cache
    t1, k1 = next(sampler, 0.0, rng)

    # Second call returns same result (idempotent)
    t2, k2 = next(sampler, 0.0, rng)
    @test t1 == t2
    @test k1 == k2
end


@safetestset rssa_bound_factor_constructor = "RSSA bound_factor in constructor" begin
    using CompetingClocks: RSSA

    # Test default bound_factor
    s1 = RSSA{Int,Float64}()
    @test s1.bound_factor == 1.05

    # Test custom bound_factor
    s2 = RSSA{Int,Float64}(bound_factor=1.2)
    @test s2.bound_factor == 1.2

    # Test bound_factor clamped to >= 1.0
    s3 = RSSA{Int,Float64}(bound_factor=0.5)
    @test s3.bound_factor == 1.0
end


@safetestset rssa_disable_errors = "RSSA disable errors" begin
    using CompetingClocks: RSSA, enable!, disable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(12345)
    sampler = RSSA{Int,Float64}()

    # Disable non-existent clock
    @test_throws KeyError disable!(sampler, 1, 0.0)

    # Enable then disable
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    disable!(sampler, 1, 0.0)

    # Disable already disabled clock
    @test_throws KeyError disable!(sampler, 1, 0.0)
end


@safetestset rssa_fire = "RSSA fire!" begin
    using CompetingClocks: RSSA, enable!, fire!, next
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(23456)
    sampler = RSSA{Int,Float64}()
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(1.0), 0.0, 0.0, rng)

    @test length(sampler) == 2

    # Get next event
    t, k = next(sampler, 0.0, rng)

    # Fire removes the clock
    fire!(sampler, k, t)
    @test length(sampler) == 1
    @test !haskey(sampler, k)
end
