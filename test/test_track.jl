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

@safetestset track_trajectory_comparison = "TrajectoryWatcher comparison" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(BigInt(2576889945234392378678934324582349))
    watcher = TrajectoryWatcher{Int64,Float64}()
    # This test will create a chain where there is one transition enabled
    # and that one transition fires, then the next is enabled.
    # At the end, it compares the trajectory likelihood with the step likelihood.
    running_loglikelihood = zero(Float64)
    check = zero(Float64)
    # Each loop represents a duration, so it creates an event and fires it.
    for step_idx in 1:10
        period_start = Float64(step_idx - 1)
        period_finish = Float64(step_idx)
        dist = Exponential(step_idx)
        λ = inv(dist.θ)
        enable!(watcher, step_idx, dist, period_start, period_start, rng)
        step_ll = steploglikelihood(watcher, period_start, period_finish, step_idx)
        @test step_ll <= 0.0
        fire!(watcher, step_idx, period_finish)
        check_step = log(λ) - λ * (period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = pathloglikelihood(watcher, Float64(10))
    @test abs(running_loglikelihood - ll) < 1e-6
end


@safetestset track_trajectory_gamma = "TrajectoryWatcher compare gamma" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(BigInt(2123459945234392378678934324582349))
    watcher = TrajectoryWatcher{Int64,Float64}()
    # As above, but this uses a non-Exponential distribution for each step.
    running_loglikelihood = zero(Float64)
    check = zero(Float64)
    # Each loop represents a duration, so it creates an event and fires it.
    for step_idx in 1:10
        period_start = Float64(step_idx - 1)
        period_finish = Float64(step_idx)
        dist = Gamma(step_idx, 1.5)
        enable!(watcher, step_idx, dist, period_start, period_start, rng)
        step_ll = steploglikelihood(watcher, period_start, period_finish, step_idx)
        @test step_ll <= 0.0
        fire!(watcher, step_idx, period_finish)
        check_step = loglikelihood(dist, period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = pathloglikelihood(watcher, Float64(10))
    @test abs(running_loglikelihood - ll) < 1e-6
end


@safetestset track_trajectory_compete = "TrajectoryWatcher compare compete" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(BigInt(2123459945234392378678934324582349))
    watcher = TrajectoryWatcher{Int64,Float64}()
    # Reset at each time step, but have one competing process in the likelihood.
    running_loglikelihood = zero(Float64)
    check = zero(Float64)
    compete = Exponential(1.0)
    compete_idx = 1003
    # Each loop represents a duration, so it creates an event and fires it.
    for step_idx in 1:10
        period_start = Float64(step_idx - 1)
        period_finish = Float64(step_idx)
        dist = Gamma(step_idx, 1.5)
        enable!(watcher, step_idx, dist, period_start, period_start, rng)
        enable!(watcher, compete_idx, compete, period_start, period_start, rng)
        step_ll = steploglikelihood(watcher, period_start, period_finish, step_idx)
        @test step_ll <= 0.0
        fire!(watcher, step_idx, period_finish)
        disable!(watcher, compete_idx, period_finish)
        check_step = loglikelihood(dist, period_finish - period_start)
        check_step += logccdf(compete, period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = pathloglikelihood(watcher, Float64(10))
    @test abs(running_loglikelihood - ll) < 1e-6
end


@safetestset track_trajectory_free = "TrajectoryWatcher compare free" begin
    using Distributions
    using DataStructures
    using CompetingClocks
    using Random
    rng = Xoshiro(BigInt(212345994379212378678934324582349))
    watcher = TrajectoryWatcher{Int64,Float64}()
    # Make a more general set of transitions and compare step-wise with total
    # log-likelihood. Make 10 transitions. At each step, fire 1, disable 4,
    # and enable 5.
    enabled = Deque{Int}()
    curtime = zero(Float64)
    for initial_idx in 1:10
        dist = Gamma(6 + initial_idx * 0.2, 1.0)
        enable!(watcher, initial_idx, dist, curtime, curtime, rng)
        push!(enabled, initial_idx)
    end
    max_idx = 10
    @test length(watcher) == max_idx

    running_loglikelihood = zero(Float64)
    for step_idx in 1:10
        tofire = popfirst!(enabled)
        firetime = curtime + 0.2
        step_ll = steploglikelihood(watcher, curtime, firetime, tofire)
        @test step_ll <= 0.0
        running_loglikelihood += step_ll
        fire!(watcher, tofire, firetime)
        curtime = firetime
        for i in 1:4
            disable_idx = popfirst!(enabled)
            disable!(watcher, disable_idx, curtime)
        end
        @test length(watcher) == 5

        # The sum of the stepwise and the trajectory-based should be equal
        # the whole way through.
        ll = pathloglikelihood(watcher, firetime)
        @test abs(running_loglikelihood - ll) < 1e-4 * abs(ll)

        for i in 1:5
            dist = Gamma(5 + 0.2 * i - 0.1 * step_idx, 1.5)
            add_idx = max_idx + 1
            enable!(watcher, add_idx, dist, curtime, curtime, rng)
            push!(enabled, add_idx)
            max_idx = add_idx
        end
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
