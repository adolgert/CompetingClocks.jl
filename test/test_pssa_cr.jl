using SafeTestsets


@safetestset pssacr_smoke = "PSSACR smoke" begin
    using CompetingClocks: PSSACR, enable!, disable!, next
    using Random: MersenneTwister
    using Distributions: Exponential

    rng = MersenneTwister(90422342)
    min_when = Inf
    seen = Set{Int}()
    sample_time = 0.5
    for i in 1:100
        sampler = PSSACR{Int,Float64}()
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


@safetestset pssacr_interface = "PSSACR basic interface" begin
    using CompetingClocks: PSSACR, enable!, disable!, next, reset!
    using Random: Xoshiro
    using Distributions: Exponential

    sampler = PSSACR{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test keytype(sampler) <: Int64

    for (clock, rate) in [(1, 0.5), (2, 1.0), (3, 2.0), (4, 10.0), (5, 0.1)]
        enable!(sampler, clock, Exponential(1/rate), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)

    disable!(sampler, 1, 0.0)
    @test !haskey(sampler, 1)
    reset!(sampler)
    @test length(sampler) == 0
end


@safetestset pssacr_empty = "PSSACR empty" begin
    using CompetingClocks: PSSACR, next
    using Random: MersenneTwister

    rng = MersenneTwister(90422342)
    sampler = PSSACR{Int,Float64}()
    when, which = next(sampler, 5.7, rng)
    @test when == Inf
    @test which === nothing
end


@safetestset pssacr_clone = "PSSACR clone" begin
    using CompetingClocks: PSSACR, enable!, clone
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(234567)
    sampler = PSSACR{Int,Float64}(ngroups=32)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)
    cloned = clone(sampler)
    @test length(cloned.groups) == 32
    @test length(cloned) == 0  # cloned is empty
    @test length(sampler) == 2  # original unchanged
end


@safetestset pssacr_copy = "PSSACR copy_clocks!" begin
    using CompetingClocks: PSSACR, enable!, next, copy_clocks!
    using Random: MersenneTwister
    using Distributions: Exponential

    src = PSSACR{Int,Float64}(ngroups=16)
    dst = PSSACR{Int,Float64}(ngroups=16)
    rng = MersenneTwister(90422342)
    enable!(src, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(src, 2, Exponential(1.0), 0.0, 0.0, rng)
    enable!(dst, 3, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(dst) == 2
    enable!(src, 5, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    enable!(dst, 6, Exponential(1.0), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3

    # Test error on mismatched ngroups
    dst_bad = PSSACR{Int,Float64}(ngroups=8)
    @test_throws ArgumentError copy_clocks!(dst_bad, src)
end


@safetestset pssacr_assign_group = "PSSACR assign_group!" begin
    using CompetingClocks: PSSACR, enable!, assign_group!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(345678)
    sampler = PSSACR{Int,Float64}(ngroups=4)

    # Assign clock 1 to group 3 before enabling
    assign_group!(sampler, 1, 3)
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    @test sampler.group_of[1] == 3

    # Test bounds checking
    @test_throws ArgumentError assign_group!(sampler, 2, 5)  # out of range
    @test_throws ArgumentError assign_group!(sampler, 2, 0)  # out of range
end


@safetestset pssacr_update_rate = "PSSACR rate update" begin
    using CompetingClocks: PSSACR, enable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = PSSACR{Int,Float64}(ngroups=1)  # single group for easy tracking

    # Enable with rate 1.0
    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0
    @test sampler.rates[1] == 1.0
    g = sampler.group_of[1]
    @test sampler.group_max[g] == 1.0

    # Update to higher rate (new_rate > group_max)
    enable!(sampler, 1, Exponential(0.5), 0.0, 0.0, rng)  # rate = 2.0
    @test sampler.rates[1] == 2.0
    @test sampler.group_max[g] == 2.0

    # Add another clock with lower rate
    enable!(sampler, 2, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0
    @test sampler.group_max[g] == 2.0  # max unchanged

    # Update clock 1 to lower rate (was max, triggers _recompute_group_max!)
    enable!(sampler, 1, Exponential(10.0), 0.0, 0.0, rng)  # rate = 0.1
    @test sampler.rates[1] == 0.1
    @test sampler.group_max[g] == 1.0  # now clock 2 is max
end


@safetestset pssacr_non_exponential = "PSSACR non-exponential error" begin
    using CompetingClocks: PSSACR, enable!
    using Random: Xoshiro
    using Distributions: Gamma, Weibull

    rng = Xoshiro(567890)
    sampler = PSSACR{Int,Float64}()

    @test_throws ArgumentError enable!(sampler, 1, Gamma(1.0, 1.0), 0.0, 0.0, rng)
    @test_throws ArgumentError enable!(sampler, 2, Weibull(1.0, 1.0), 0.0, 0.0, rng)
end


@safetestset pssacr_swap_remove = "PSSACR swap-remove in disable" begin
    using CompetingClocks: PSSACR, enable!, disable!
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(678901)
    # Use single group to ensure all clocks go to same group
    sampler = PSSACR{Int,Float64}(ngroups=1)

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 3, Exponential(1.0), 0.0, 0.0, rng)

    # Clock 1 is first in group, clock 3 is last
    # Disabling clock 1 should trigger swap-remove (pos != lastidx)
    disable!(sampler, 1, 0.0)

    @test !haskey(sampler, 1)
    @test haskey(sampler, 2)
    @test haskey(sampler, 3)
    @test length(sampler) == 2
end


@safetestset pssacr_getindex = "PSSACR getindex" begin
    using CompetingClocks: PSSACR, enable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    rng = Xoshiro(789012)
    sampler = PSSACR{Int,Float64}()
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
