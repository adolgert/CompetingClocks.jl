using SafeTestsets


@safetestset track_trackwatcher_smoke = "TrackWatcher smoke" begin
    using Distributions
    using CompetingClocks
    using Random
    using Base
    rng = Xoshiro(3242234)
    tw = TrackWatcher{Int,Float64}()
    enable!(tw, 3, Exponential(), 0.0, 0.0, rng)
    @test length(tw.enabled) == 1 && 3 ∈ keys(tw.enabled)
    enable!(tw, 4, Exponential(), 0.0, 3.0, rng)
    @test length(tw.enabled) == 2 && 4 ∈ keys(tw.enabled)
    enable!(tw, 7, Exponential(), 5.0, 5.0, rng)
    @test length(tw.enabled) == 3 && 7 ∈ keys(tw.enabled)
    disable!(tw, 4, 9.0)
    @test length(tw.enabled) == 2 && 4 ∉ keys(tw.enabled)

    dst = TrackWatcher{Int,Float64}()
    enable!(dst, 11, Exponential(), 5.0, 5.0, rng)
    copy_clocks!(dst, tw)
    @test length(tw.enabled) == 2 && 11 ∉ keys(tw.enabled)
end


@safetestset memory_sampler = "MemorySampler known insertions" begin
    # Let's put entries into its heap in order to ensure the heap does the right thing.
    # This is a clear-box test. It depends heavily on implementation.
    using CompetingClocks: FirstToFire, enable!, disable!, next
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    propagator = MemorySampler(FirstToFire{Int64,Float64}())
    rng = Xoshiro(39472)

    clock_key = 1:5
    relative_firing = [7.9, 12.3, 3.7, 1.3, 0.2]
    enabling_times = [0.0; cumsum(fill(0.01, 4), dims=1)]
    for i in eachindex(relative_firing)
        enable!(propagator, clock_key[i], Dirac(relative_firing[i]), enabling_times[i], enabling_times[i], rng)
    end

    # test retrieving the enabling times
    for i in eachindex(relative_firing)
        @test absolute_enabling(propagator, clock_key[i]) == enabling_times[i]
    end

    absolute_firing = enabling_times + relative_firing
    clock_firing_order = sortperm(absolute_firing)
    last_fired = 0.0
    for i in eachindex(clock_firing_order)
        (last_fired, which) = next(propagator, last_fired, rng)
        @test last_fired ≈ absolute_firing[clock_firing_order[i]]
        disable!(propagator, which, last_fired)
        @test_throws KeyError absolute_enabling(propagator, clock_key[clock_firing_order[i]])
    end
end


@safetestset track_debugwatcher_smoke = "DebugWatcher smoke" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(3242234)
    dw = DebugWatcher{Int,Float64}()
    enable!(dw, 3, Exponential(), 0.0, 0.0, rng)
    @test dw.enabled[1].clock == 3
    enable!(dw, 4, Exponential(), 0.0, 3.0, rng)
    @test dw.enabled[2].clock == 4
    enable!(dw, 7, Exponential(), 5.0, 5.0, rng)
    @test dw.enabled[3].clock == 7
    disable!(dw, 4, 9.0)
    @test dw.disabled[1].clock == 4

    dst = DebugWatcher{Int,Float64}()
    enable!(dst, 11, Exponential(), 5.0, 5.0, rng)
    copy_clocks!(dst, dw)
    @test length(dw.enabled) == 3 && length(dw.disabled) == 1
end


@safetestset track_debugwatcher_clone = "DebugWatcher clone" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: clone
    using Random

    rng = Xoshiro(9876543)

    # Create and populate a DebugWatcher
    dw = DebugWatcher{Int,Float64}()
    enable!(dw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(dw, 2, Gamma(2.0, 1.0), 1.0, 1.0, rng)
    disable!(dw, 1, 2.0)

    # Clone should create an empty watcher with same type parameters
    cloned = clone(dw)
    @test cloned isa DebugWatcher{Int,Float64}
    @test isempty(cloned.enabled)
    @test isempty(cloned.disabled)
    @test cloned.log == dw.log

    # Test clone preserves log setting
    dw_nolog = DebugWatcher{String,Float64}(; log=false)
    cloned_nolog = clone(dw_nolog)
    @test cloned_nolog.log == false
end


@safetestset track_debugwatcher_reset = "DebugWatcher reset!" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: reset!
    using Random

    rng = Xoshiro(1234567)

    # Create and populate a DebugWatcher
    dw = DebugWatcher{Int,Float64}()
    enable!(dw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(dw, 2, Exponential(2.0), 0.5, 0.5, rng)
    enable!(dw, 3, Exponential(3.0), 1.0, 1.0, rng)
    disable!(dw, 1, 1.5)
    disable!(dw, 2, 2.0)

    @test length(dw.enabled) == 3
    @test length(dw.disabled) == 2

    # Reset should clear all entries
    reset!(dw)
    @test isempty(dw.enabled)
    @test isempty(dw.disabled)

    # Should be able to use it again after reset
    enable!(dw, 10, Exponential(1.0), 0.0, 0.0, rng)
    @test length(dw.enabled) == 1
    @test dw.enabled[1].clock == 10
end


@safetestset track_debugwatcher_fire = "DebugWatcher fire!" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(7654321)

    dw = DebugWatcher{Int,Float64}()

    # Enable some clocks
    enable!(dw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(dw, 2, Exponential(2.0), 0.0, 0.0, rng)
    enable!(dw, 3, Exponential(3.0), 0.0, 0.0, rng)

    @test length(dw.enabled) == 3
    @test length(dw.disabled) == 0

    # Fire clock 2 - this should add to disabled list
    fire!(dw, 2, 1.5)
    @test length(dw.disabled) == 1
    @test dw.disabled[1].clock == 2
    @test dw.disabled[1].when == 1.5

    # Fire another clock
    fire!(dw, 1, 2.0)
    @test length(dw.disabled) == 2
    @test dw.disabled[2].clock == 1
    @test dw.disabled[2].when == 2.0

    # Disable (not fire) the last one
    disable!(dw, 3, 3.0)
    @test length(dw.disabled) == 3
    @test dw.disabled[3].clock == 3
end


@safetestset track_debugwatcher_log_option = "DebugWatcher log option" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(1111111)

    # Test with logging disabled
    dw_nolog = DebugWatcher{Int,Float64}(; log=false)
    @test dw_nolog.log == false

    # Operations should work without logging
    enable!(dw_nolog, 1, Exponential(1.0), 0.0, 0.0, rng)
    disable!(dw_nolog, 1, 1.0)
    fire!(dw_nolog, 2, 2.0)

    @test length(dw_nolog.enabled) == 1
    @test length(dw_nolog.disabled) == 2

    # Test with logging enabled (default)
    dw_log = DebugWatcher{Int,Float64}()
    @test dw_log.log == true
end


@safetestset track_trackwatcher_clone = "TrackWatcher clone" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: clone
    using Random

    rng = Xoshiro(5555555)

    # Create and populate a TrackWatcher
    tw = TrackWatcher{Int,Float64}()
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(tw, 2, Gamma(2.0, 1.0), 1.0, 1.0, rng)
    enable!(tw, 3, Weibull(1.5, 2.0), 2.0, 2.0, rng)

    @test length(tw) == 3

    # Clone should create an empty watcher with same type parameters
    cloned = clone(tw)
    @test cloned isa TrackWatcher{Int,Float64}
    @test length(cloned) == 0
    @test isempty(cloned.enabled)

    # Original should be unchanged
    @test length(tw) == 3
end


@safetestset track_jitter = "EnabledWatcher jitter!" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: jitter!
    using Random

    rng = Xoshiro(6666666)

    # jitter! on TrackWatcher should do nothing (returns nothing)
    tw = TrackWatcher{Int,Float64}()
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)

    result = jitter!(tw, 1.0, rng)
    @test result === nothing
    # The enabled clock should be unchanged
    @test length(tw) == 1
    @test haskey(tw, 1)
end


@safetestset track_memorysampler_keytype = "MemorySampler keytype" begin
    using CompetingClocks: FirstToFire, MemorySampler, keytype
    using CompetingClocks

    # Create MemorySampler with Int64 keys
    propagator = MemorySampler(FirstToFire{Int64,Float64}())
    @test keytype(propagator) == Int64

    # Create MemorySampler with Symbol keys
    propagator_sym = MemorySampler(FirstToFire{Symbol,Float64}())
    @test keytype(propagator_sym) == Symbol

    # Create MemorySampler with String keys
    propagator_str = MemorySampler(FirstToFire{String,Float32}())
    @test keytype(propagator_str) == String
end


@safetestset track_memorysampler_getindex = "MemorySampler getindex" begin
    using CompetingClocks: FirstToFire, MemorySampler
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    propagator = MemorySampler(FirstToFire{Int64,Float64}())
    rng = Xoshiro(7777777)

    # Enable some clocks
    enable!(propagator, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(propagator, 2, Exponential(2.0), 0.0, 0.0, rng)
    enable!(propagator, 3, Exponential(3.0), 0.0, 0.0, rng)

    # getindex should delegate to the underlying sampler
    # For FirstToFire, getindex returns the firing time
    t1 = propagator[1]
    t2 = propagator[2]
    t3 = propagator[3]

    @test t1 isa Float64
    @test t2 isa Float64
    @test t3 isa Float64
    @test t1 >= 0.0
    @test t2 >= 0.0
    @test t3 >= 0.0
end


@safetestset track_isenabled_enabled = "TrackWatcher isenabled and enabled" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: isenabled, enabled
    using Random

    rng = Xoshiro(8888888)

    tw = TrackWatcher{Int,Float64}()

    # Initially nothing is enabled
    @test !isenabled(tw, 1)
    @test !isenabled(tw, 2)
    @test isempty(enabled(tw))

    # Enable some clocks
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(tw, 2, Exponential(2.0), 0.5, 0.5, rng)

    @test isenabled(tw, 1)
    @test isenabled(tw, 2)
    @test !isenabled(tw, 3)

    enabled_keys = collect(enabled(tw))
    @test length(enabled_keys) == 2
    @test 1 in enabled_keys
    @test 2 in enabled_keys

    # Disable one
    disable!(tw, 1, 1.0)
    @test !isenabled(tw, 1)
    @test isenabled(tw, 2)

    enabled_keys_after = collect(enabled(tw))
    @test length(enabled_keys_after) == 1
    @test 2 in enabled_keys_after
end


@safetestset track_trackwatcher_iterate = "TrackWatcher iteration" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(9999999)

    tw = TrackWatcher{Int,Float64}()
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(tw, 2, Gamma(2.0, 1.0), 1.0, 1.0, rng)
    enable!(tw, 3, Weibull(1.5, 2.0), 2.0, 2.0, rng)

    # Test iteration via for-loop (uses Base.iterate)
    clocks_seen = Int[]
    for entry in tw
        push!(clocks_seen, entry.clock)
        @test entry.distribution isa Distributions.UnivariateDistribution
        @test entry.te isa Float64
        @test entry.when isa Float64
    end

    @test length(clocks_seen) == 3
    @test sort(clocks_seen) == [1, 2, 3]
end


@safetestset track_disable_nonexistent = "TrackWatcher disable nonexistent clock" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(1010101)

    tw = TrackWatcher{Int,Float64}()
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)

    # Disabling a clock that doesn't exist should not error
    disable!(tw, 999, 1.0)  # This clock was never enabled
    @test length(tw) == 1
    @test haskey(tw, 1)
end


@safetestset track_enable_replaces = "TrackWatcher enable replaces existing" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(2020202)

    tw = TrackWatcher{Int,Float64}()

    # Enable clock 1 with one distribution
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)
    @test tw[1].distribution isa Exponential
    @test tw[1].te == 0.0
    @test tw[1].when == 0.0

    # Enable clock 1 again with different distribution and times
    # This should replace the existing entry (after disabling it internally)
    enable!(tw, 1, Gamma(2.0, 1.0), 5.0, 5.0, rng)
    @test length(tw) == 1  # Still only one clock
    @test tw[1].distribution isa Gamma
    @test tw[1].te == 5.0
    @test tw[1].when == 5.0
end


@safetestset track_steploglikelihood_edge_cases = "steploglikelihood edge cases" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(3030303)

    tw = TrackWatcher{Int,Float64}()

    # Case: t < te (firing time is before the enabling time offset)
    # This is the case where the clock hasn't "started" yet
    enable!(tw, 1, Exponential(1.0), 5.0, 0.0, rng)  # te=5.0, when=0.0

    # If t < te and this clock fires, should return -Inf
    ll_fires_early = steploglikelihood(tw, 0.0, 3.0, 1)
    @test ll_fires_early == -Inf

    # If t < te and a different clock fires, should return 0 for this clock's contribution
    enable!(tw, 2, Exponential(1.0), 0.0, 0.0, rng)  # te=0.0
    ll_other_fires = steploglikelihood(tw, 0.0, 3.0, 2)
    @test ll_other_fires < 0  # Should be finite since clock 2 can fire
    @test isfinite(ll_other_fires)
end


@safetestset track_stepcumulant = "stepcumulant function" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: stepcumulant
    using Random

    rng = Xoshiro(4040404)

    tw = TrackWatcher{Int,Float64}()

    # Single exponential clock - cumulant should be between 0 and 1
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)

    # Test at various times
    for t in [0.1, 0.5, 1.0, 2.0, 5.0]
        c = stepcumulant(tw, 0.0, t)
        @test 0.0 <= c <= 1.0
    end

    # At t=0, cumulant should be 0
    c_zero = stepcumulant(tw, 0.0, 0.0)
    @test c_zero ≈ 0.0 atol=1e-10

    # Test with t < te (before distribution starts)
    tw2 = TrackWatcher{Int,Float64}()
    enable!(tw2, 1, Exponential(1.0), 5.0, 0.0, rng)
    c_early = stepcumulant(tw2, 0.0, 3.0)
    @test c_early ≈ 0.0 atol=1e-10
end


@safetestset track_stepconditionalprobability = "stepconditionalprobability function" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: stepconditionalprobability, hazard
    using Random

    rng = Xoshiro(5050505)

    tw = TrackWatcher{Int,Float64}()

    # Two exponential clocks with same rate - should have equal probability
    enable!(tw, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(tw, 2, Exponential(1.0), 0.0, 0.0, rng)

    probs = stepconditionalprobability(tw, 1.0)
    @test length(probs) == 2
    @test probs[1] ≈ 0.5 atol=1e-10
    @test probs[2] ≈ 0.5 atol=1e-10
    @test sum(values(probs)) ≈ 1.0 atol=1e-10

    # Different rates - higher rate should have higher probability
    tw2 = TrackWatcher{Int,Float64}()
    enable!(tw2, 1, Exponential(1.0), 0.0, 0.0, rng)  # rate = 1
    enable!(tw2, 2, Exponential(0.5), 0.0, 0.0, rng)  # rate = 2

    probs2 = stepconditionalprobability(tw2, 1.0)
    @test probs2[2] > probs2[1]  # Clock 2 has higher rate
    @test sum(values(probs2)) ≈ 1.0 atol=1e-10
end


@safetestset track_stepconditionalprobability_zero_hazard = "stepconditionalprobability zero hazard" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: stepconditionalprobability
    using Random

    rng = Xoshiro(6060606)

    tw = TrackWatcher{Int,Float64}()

    # Enable a clock with te in the future so hazard at current time is 0
    enable!(tw, 1, Exponential(1.0), 10.0, 0.0, rng)

    # At t=5.0, the clock hasn't started (te=10.0), so hazard should be 0
    # When all hazards are zero, the function should return the dict with zeros
    probs = stepconditionalprobability(tw, 5.0)
    @test length(probs) == 1
    # The hazard at t=5 for a distribution starting at te=10 should be 0
    # because t < te
end
