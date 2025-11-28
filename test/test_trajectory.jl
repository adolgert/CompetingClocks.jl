using SafeTestsets

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


@safetestset trajectory_clone = "TrajectoryWatcher clone" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(1234)
    watcher = TrajectoryWatcher{Int64,Float64}()

    # Add some state
    enable!(watcher, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(watcher, 2, Gamma(2.0, 1.5), 0.0, 0.0, rng)
    fire!(watcher, 1, 0.5)

    # Clone should create a fresh instance
    cloned = CompetingClocks.clone(watcher)
    @test length(cloned) == 0
    @test typeof(cloned) == typeof(watcher)
end


@safetestset trajectory_copy_clocks = "TrajectoryWatcher copy_clocks!" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(5678)

    src = TrajectoryWatcher{Int64,Float64}()
    enable!(src, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(src, 2, Gamma(2.0, 1.5), 0.1, 0.1, rng)
    fire!(src, 1, 0.5)

    dst = TrajectoryWatcher{Int64,Float64}()
    enable!(dst, 99, Exponential(2.0), 0.0, 0.0, rng)

    CompetingClocks.copy_clocks!(dst, src)

    # dst should have same enabled clocks as src
    @test length(dst) == length(src)
    @test haskey(dst.enabled, 2)
    @test !haskey(dst.enabled, 99)
    @test dst.loglikelihood == src.loglikelihood
    @test dst.curtime == src.curtime
end


@safetestset trajectory_reset = "TrajectoryWatcher reset!" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(9012)

    watcher = TrajectoryWatcher{Int64,Float64}()
    enable!(watcher, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(watcher, 2, Gamma(2.0, 1.5), 0.0, 0.0, rng)
    fire!(watcher, 1, 0.5)

    @test length(watcher) > 0
    @test watcher.loglikelihood != 0.0
    @test watcher.curtime != 0.0

    CompetingClocks.reset!(watcher)

    @test length(watcher) == 0
    @test watcher.loglikelihood == 0.0
    @test watcher.curtime == 0.0
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


# ---------------------------------------------------------------------------- #
#                            PathLikelihoods tests                             #
# ---------------------------------------------------------------------------- #

@safetestset pathlikelihoods_basic = "PathLikelihoods basic operations" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(1234)

    pl = PathLikelihoods{Int64,Float64}(3)

    # Test initial state
    @test length(pl) == 0
    @test pl.loglikelihood == zeros(Float64, 3)
    @test pl.curtime == 0.0

    # Test enable! with single distribution
    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)
    @test length(pl) == 1
    @test CompetingClocks.isenabled(pl, 1)
    @test 1 in CompetingClocks.enabled(pl)

    # Test enable! overwrites existing clock
    enable!(pl, 1, Exponential(2.0), 0.1, 0.1, rng)
    @test length(pl) == 1

    # Add another clock
    enable!(pl, 2, Gamma(2.0, 1.5), 0.0, 0.0, rng)
    @test length(pl) == 2
end


@safetestset pathlikelihoods_clone = "PathLikelihoods clone" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(5678)

    pl = PathLikelihoods{Int64,Float64}(4)
    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)
    fire!(pl, 1, 0.5)

    cloned = CompetingClocks.clone(pl)
    @test length(cloned) == 0
    @test length(cloned.loglikelihood) == 4
    @test cloned.loglikelihood == zeros(Float64, 4)
end


@safetestset pathlikelihoods_copy_clocks = "PathLikelihoods copy_clocks!" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(9012)

    src = PathLikelihoods{Int64,Float64}(3)
    enable!(src, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(src, 2, Gamma(2.0, 1.5), 0.1, 0.1, rng)
    fire!(src, 1, 0.5)

    dst = PathLikelihoods{Int64,Float64}(3)
    enable!(dst, 99, Exponential(2.0), 0.0, 0.0, rng)

    CompetingClocks.copy_clocks!(dst, src)

    @test length(dst) == length(src)
    @test CompetingClocks.isenabled(dst, 2)
    @test !CompetingClocks.isenabled(dst, 99)
    @test dst.loglikelihood == src.loglikelihood
    @test dst.curtime == src.curtime
end


@safetestset pathlikelihoods_reset = "PathLikelihoods reset!" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(3456)

    pl = PathLikelihoods{Int64,Float64}(3)
    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(pl, 2, Gamma(2.0, 1.5), 0.0, 0.0, rng)
    fire!(pl, 1, 0.5)

    @test length(pl) > 0
    @test any(x -> x != 0.0, pl.loglikelihood)
    @test pl.curtime != 0.0

    CompetingClocks.reset!(pl)

    @test length(pl) == 0
    @test all(x -> x == 0.0, pl.loglikelihood)
    @test pl.curtime == 0.0
end


@safetestset pathlikelihoods_vs_trajectory_single = "PathLikelihoods matches TrajectoryWatcher with single dist" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(7890)

    # Both should give same result when PathLikelihoods uses single distributions
    tw = TrajectoryWatcher{Int64,Float64}()
    pl = PathLikelihoods{Int64,Float64}(1)

    # Run the same sequence of operations
    enable!(tw, 1, Exponential(2.0), 0.0, 0.0, rng)
    enable!(pl, 1, Exponential(2.0), 0.0, 0.0, rng)

    enable!(tw, 2, Gamma(3.0, 1.0), 0.0, 0.0, rng)
    enable!(pl, 2, Gamma(3.0, 1.0), 0.0, 0.0, rng)

    # Fire clock 1
    fire!(tw, 1, 0.5)
    fire!(pl, 1, 0.5)

    # Disable clock 2
    disable!(tw, 2, 0.5)
    disable!(pl, 2, 0.5)

    # Enable new clocks
    enable!(tw, 3, Weibull(2.0, 1.5), 0.5, 0.5, rng)
    enable!(pl, 3, Weibull(2.0, 1.5), 0.5, 0.5, rng)

    # Fire clock 3
    fire!(tw, 3, 1.0)
    fire!(pl, 3, 1.0)

    # Compare likelihoods
    tw_ll = pathloglikelihood(tw, 1.0)
    pl_ll = pathloglikelihood(pl, 1.0)

    @test length(pl_ll) == 1
    @test abs(tw_ll - pl_ll[1]) < 1e-10
end


@safetestset pathlikelihoods_vs_trajectory_chain = "PathLikelihoods matches TrajectoryWatcher chain" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(BigInt(2576889945234392378678934324582349))
    tw = TrajectoryWatcher{Int64,Float64}()
    pl = PathLikelihoods{Int64,Float64}(1)

    # Run the same chain test as track_trajectory_comparison
    for step_idx in 1:10
        period_start = Float64(step_idx - 1)
        period_finish = Float64(step_idx)
        dist = Exponential(step_idx)
        enable!(tw, step_idx, dist, period_start, period_start, rng)
        enable!(pl, step_idx, dist, period_start, period_start, rng)
        fire!(tw, step_idx, period_finish)
        fire!(pl, step_idx, period_finish)
    end

    tw_ll = pathloglikelihood(tw, Float64(10))
    pl_ll = pathloglikelihood(pl, Float64(10))

    @test abs(tw_ll - pl_ll[1]) < 1e-10
end


@safetestset pathlikelihoods_vs_trajectory_compete = "PathLikelihoods matches TrajectoryWatcher compete" begin
    using Distributions
    using CompetingClocks
    using Random

    rng = Xoshiro(BigInt(2123459945234392378678934324582349))
    tw = TrajectoryWatcher{Int64,Float64}()
    pl = PathLikelihoods{Int64,Float64}(1)

    compete = Exponential(1.0)
    compete_idx = 1003

    for step_idx in 1:10
        period_start = Float64(step_idx - 1)
        period_finish = Float64(step_idx)
        dist = Gamma(step_idx, 1.5)
        enable!(tw, step_idx, dist, period_start, period_start, rng)
        enable!(pl, step_idx, dist, period_start, period_start, rng)
        enable!(tw, compete_idx, compete, period_start, period_start, rng)
        enable!(pl, compete_idx, compete, period_start, period_start, rng)
        fire!(tw, step_idx, period_finish)
        fire!(pl, step_idx, period_finish)
        disable!(tw, compete_idx, period_finish)
        disable!(pl, compete_idx, period_finish)
    end

    tw_ll = pathloglikelihood(tw, Float64(10))
    pl_ll = pathloglikelihood(pl, Float64(10))

    @test abs(tw_ll - pl_ll[1]) < 1e-10
end


@safetestset pathlikelihoods_vector_distributions = "PathLikelihoods with vector distributions" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(1357)

    # Create PathLikelihoods with space for 3 likelihood calculations
    pl = PathLikelihoods{Int64,Float64}(3)

    # Enable with a vector of distributions (different parameters)
    dists = [Exponential(1.0), Exponential(2.0), Exponential(3.0)]
    enable!(pl, 1, dists, 0.0, 0.0, rng)

    @test length(pl) == 1

    # Fire the clock
    fire!(pl, 1, 0.5)

    # Each likelihood should be different because the distributions differ
    ll = pathloglikelihood(pl, 0.5)
    @test length(ll) == 3

    # Manually check the likelihoods: logpdf(Exponential(λ), t) = log(1/λ) - t/λ
    for idx in 1:3
        λ = Float64(idx)  # 1.0, 2.0, 3.0
        expected = log(1/λ) - 0.5/λ
        @test abs(ll[idx] - expected) < 1e-10
    end

    # Verify they're all different
    @test ll[1] != ll[2]
    @test ll[2] != ll[3]
end


@safetestset pathlikelihoods_vector_disable = "PathLikelihoods vector distributions disable" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(2468)

    pl = PathLikelihoods{Int64,Float64}(3)

    # Enable with vector of distributions
    dists = [Exponential(1.0), Exponential(2.0), Exponential(3.0)]
    enable!(pl, 1, dists, 0.0, 0.0, rng)

    # Disable without firing
    disable!(pl, 1, 0.5)

    # Likelihood should include the survival (ccdf) for each distribution
    ll = pathloglikelihood(pl, 0.5)
    @test length(ll) == 3

    # logccdf(Exponential(λ), t) = -t/λ
    for idx in 1:3
        λ = Float64(idx)
        expected = -0.5/λ
        @test abs(ll[idx] - expected) < 1e-10
    end
end


@safetestset pathlikelihoods_mixed_single_vector = "PathLikelihoods mixed single and vector" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(3579)

    pl = PathLikelihoods{Int64,Float64}(3)

    # Enable clock 1 with single distribution (should broadcast to all 3)
    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)

    # Enable clock 2 with vector of distributions
    dists = [Gamma(2.0, 1.0), Gamma(2.0, 1.5), Gamma(2.0, 2.0)]
    enable!(pl, 2, dists, 0.0, 0.0, rng)

    @test length(pl) == 2

    # Fire clock 1
    fire!(pl, 1, 0.3)

    # Disable clock 2
    disable!(pl, 2, 0.3)

    ll = pathloglikelihood(pl, 0.3)
    @test length(ll) == 3

    # All should include the same log-pdf for clock 1
    exp_ll_clock1 = logpdf(Exponential(1.0), 0.3)

    # Each should have different ccdf for clock 2
    for idx in 1:3
        scale = [1.0, 1.5, 2.0][idx]
        exp_ll_clock2 = logccdf(Gamma(2.0, scale), 0.3)
        expected = exp_ll_clock1 + exp_ll_clock2
        @test abs(ll[idx] - expected) < 1e-10
    end
end


@safetestset pathlikelihoods_shifted_te = "PathLikelihoods with shifted enabling time" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(4680)

    # Test the case where te < when (distribution zero-point is before enabling time)
    tw = TrajectoryWatcher{Int64,Float64}()
    pl = PathLikelihoods{Int64,Float64}(1)

    # Enable with te=0.0 but when=0.5 (distribution started earlier)
    dist = Exponential(2.0)
    enable!(tw, 1, dist, 0.0, 0.5, rng)
    enable!(pl, 1, dist, 0.0, 0.5, rng)

    # Fire at time 1.0
    fire!(tw, 1, 1.0)
    fire!(pl, 1, 1.0)

    tw_ll = pathloglikelihood(tw, 1.0)
    pl_ll = pathloglikelihood(pl, 1.0)

    @test abs(tw_ll - pl_ll[1]) < 1e-10

    # The likelihood should account for the shifted enabling:
    # logpdf at (1.0 - 0.0) = 1.0, minus logccdf at (0.5 - 0.0) = 0.5
    expected = logpdf(dist, 1.0) - logccdf(dist, 0.5)
    @test abs(tw_ll - expected) < 1e-10
end


@safetestset pathlikelihoods_remaining_clocks = "PathLikelihoods remaining clocks in pathloglikelihood" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(5791)

    # Test that pathloglikelihood properly handles clocks that haven't fired yet
    pl = PathLikelihoods{Int64,Float64}(1)

    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(pl, 2, Gamma(2.0, 1.5), 0.0, 0.0, rng)

    # Fire only clock 1, leave clock 2 enabled
    fire!(pl, 1, 0.5)

    # pathloglikelihood at time 1.0 should include survival of clock 2
    ll = pathloglikelihood(pl, 1.0)

    # Expected: logpdf(Exp(1), 0.5) + logccdf(Gamma(2,1.5), 1.0)
    expected = logpdf(Exponential(1.0), 0.5) + logccdf(Gamma(2.0, 1.5), 1.0)
    @test abs(ll[1] - expected) < 1e-10
end


@safetestset pathlikelihoods_pathloglikelihood_before_curtime = "PathLikelihoods pathloglikelihood when <= curtime" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(6802)

    pl = PathLikelihoods{Int64,Float64}(1)

    enable!(pl, 1, Exponential(1.0), 0.0, 0.0, rng)
    fire!(pl, 1, 0.5)

    # curtime is now 0.5
    # Calling pathloglikelihood with when <= curtime should just return loglikelihood
    ll_at_curtime = pathloglikelihood(pl, 0.5)
    ll_before = pathloglikelihood(pl, 0.3)

    @test ll_at_curtime == pl.loglikelihood
    @test ll_before == pl.loglikelihood
end


@safetestset pathlikelihoods_error_disable_not_enabled = "PathLikelihoods error on disable not enabled" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(7913)

    pl = PathLikelihoods{Int64,Float64}(1)

    # Try to disable a clock that was never enabled
    @test_throws ErrorException disable!(pl, 999, 0.5)
end


@safetestset pathlikelihoods_error_fire_not_enabled = "PathLikelihoods error on fire not enabled" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(8024)

    pl = PathLikelihoods{Int64,Float64}(1)

    # Try to fire a clock that was never enabled
    @test_throws ErrorException fire!(pl, 999, 0.5)
end


@safetestset trajectory_error_disable_not_enabled = "TrajectoryWatcher error on disable not enabled" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(9135)

    tw = TrajectoryWatcher{Int64,Float64}()

    # Try to disable a clock that was never enabled
    @test_throws ErrorException disable!(tw, 999, 0.5)
end


@safetestset trajectory_error_fire_not_enabled = "TrajectoryWatcher error on fire not enabled" begin
    using Distributions
    using CompetingClocks
    using Random
    rng = Xoshiro(1024)

    tw = TrajectoryWatcher{Int64,Float64}()

    # Try to fire a clock that was never enabled
    @test_throws ErrorException fire!(tw, 999, 0.5)
end
