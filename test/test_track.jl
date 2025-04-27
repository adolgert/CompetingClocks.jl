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
    copy!(dst, tw)
    @test length(tw.enabled) == 2 && 11 ∉ keys(tw.enabled)
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
        fire!(watcher, step_idx, period_finish)
        check_step = log(λ) - λ * (period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = trajectoryloglikelihood(watcher)
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
        fire!(watcher, step_idx, period_finish)
        check_step = loglikelihood(dist, period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = trajectoryloglikelihood(watcher)
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
        fire!(watcher, step_idx, period_finish)
        disable!(watcher, compete_idx, period_finish)
        check_step = loglikelihood(dist, period_finish - period_start)
        check_step += logccdf(compete, period_finish - period_start)
        @test abs(check_step - step_ll) < 1e-4 * abs(step_ll)
        running_loglikelihood += step_ll
        check += check_step
    end
    ll = trajectoryloglikelihood(watcher)
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
        ll = trajectoryloglikelihood(watcher)
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
    copy!(dst, dw)
    @test length(dw.enabled) == 3 && length(dw.disabled) == 1
end
