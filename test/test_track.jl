using SafeTestsets


@safetestset track_trackwatcher_smoke = "TrackWatcher smoke" begin
    using Distributions
    using CompetingClocks
    using Random
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
    running = zero(Float64)
    for step_idx in 1:10
        hazard = 1.0 / step_idx
        dist = Exponential(step_idx)
        enable!(watcher, step_idx, dist, Float64(step_idx - 1), Float64(step_idx - 1), rng)
        running += steploglikelihood(watcher, Float64(step_idx - 1), Float64(step_idx), step_idx)
        fire!(watcher, step_idx, Float64(step_idx))
        disable!(watcher, step_idx, Float64(step_idx))
    end
    ll = trajectoryloglikelihood(watcher)
    @test abs(running - ll) < 1e-6
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
