using SafeTestsets

# These tests pin the invariant that ForwardDiff.Dual parameters flow through
# the log-likelihood accumulators, so that gradients of a path log-likelihood
# with respect to distribution parameters can be taken by automatic
# differentiation. The watchers never draw from the rng, so a fixed seed keeps
# the traces deterministic; the rng argument is only there for interface parity.

@safetestset autodiff_exponential_gradient_matches_analytic = "TrajectoryWatcher gradient of exponential race equals analytic score" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: TrajectoryWatcher
    using ForwardDiff
    using Random
    rng = Xoshiro(1)

    # A fixed race between two exponential clocks. Each clock is re-enabled fresh
    # the instant it fires, so both are effectively enabled over the whole
    # interval [0, t_N]. Distributions.Exponential takes the SCALE = 1/rate.
    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    function loglik(θ)
        la, lb = θ[1], θ[2]
        tw = TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0, rng)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0, rng)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t, rng)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    θ = [0.7, 0.4]
    g = ForwardDiff.gradient(loglik, θ)

    # For competing exponentials always enabled over [0, t_N], the score is
    # n_k/λ_k - t_N: three :a fires, two :b fires, final time 4.5.
    analytic = [3 / 0.7 - 4.5, 2 / 0.4 - 4.5]
    @test isapprox(g, analytic; rtol=1e-10, atol=1e-12)
end


@safetestset autodiff_trajectory_matches_stepwise_gradient = "Accumulating and stepwise likelihood gradients agree under AD" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: TrajectoryWatcher, TrackWatcher
    using ForwardDiff
    using Random
    rng = Xoshiro(1)

    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    # Path A: the accumulating TrajectoryWatcher.
    function loglik_accum(θ)
        la, lb = θ[1], θ[2]
        tw = TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0, rng)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0, rng)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t, rng)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    # Path B: sum steploglikelihood over a bare TrackWatcher on the same trace.
    function loglik_step(θ)
        la, lb = θ[1], θ[2]
        track = TrackWatcher{Symbol,Float64}()
        enable!(track, :a, Exponential(1 / la), 0.0, 0.0, rng)
        enable!(track, :b, Exponential(1 / lb), 0.0, 0.0, rng)
        total = zero(eltype(θ))
        t0 = 0.0
        for (t, k) in trace
            total += steploglikelihood(track, t0, t, k)
            fire!(track, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(track, k, d, t, t, rng)
            t0 = t
        end
        return total
    end

    θ = [0.7, 0.4]
    @test isapprox(loglik_accum(θ), loglik_step(θ); rtol=1e-12, atol=1e-12)
    @test isapprox(ForwardDiff.gradient(loglik_accum, θ),
                   ForwardDiff.gradient(loglik_step, θ); rtol=1e-10, atol=1e-12)
end


@safetestset autodiff_weibull_gradient_matches_finite_difference = "TrajectoryWatcher gradient of Weibull race matches finite differences" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: TrajectoryWatcher
    using ForwardDiff
    using Random
    rng = Xoshiro(1)

    # A non-exponential clock (Weibull, shape and scale) racing an exponential.
    trace = [(0.5, :e), (1.2, :w), (2.0, :e)]

    function loglik(θ)
        shape, scale, rate = θ[1], θ[2], θ[3]
        tw = TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
        enable!(tw, :w, Weibull(shape, scale), 0.0, 0.0, rng)
        enable!(tw, :e, Exponential(1 / rate), 0.0, 0.0, rng)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :w ? Weibull(shape, scale) : Exponential(1 / rate)
            enable!(tw, k, d, t, t, rng)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    θ = [1.5, 2.0, 0.8]
    g_ad = ForwardDiff.gradient(loglik, θ)

    h = 1e-6
    g_fd = similar(θ)
    for i in eachindex(θ)
        θp = copy(θ); θp[i] += h
        θm = copy(θ); θm[i] -= h
        g_fd[i] = (loglik(θp) - loglik(θm)) / (2h)
    end

    @test isapprox(g_ad, g_fd; rtol=1e-6, atol=1e-8)
end


@safetestset autodiff_pathlikelihoods_component_matches_trajectory = "PathLikelihoods component derivative equals TrajectoryWatcher derivative" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: PathLikelihoods, TrajectoryWatcher
    using ForwardDiff
    using Random
    rng = Xoshiro(1)

    # PathLikelihoods carries a vector of candidate distributions; the first
    # component uses rate θ and must reproduce a TrajectoryWatcher on that same
    # distribution, both in value and derivative.
    function pl_component(θ, idx)
        pl = PathLikelihoods{Symbol,Float64,typeof(θ)}(2)
        dists = [Exponential(1 / θ), Exponential(1 / (2θ))]
        enable!(pl, :a, dists, 0.0, 0.0, rng)
        fire!(pl, :a, 0.7)
        return pathloglikelihood(pl, 0.7)[idx]
    end

    function tw_value(θ)
        tw = TrajectoryWatcher{Symbol,Float64,typeof(θ)}()
        enable!(tw, :a, Exponential(1 / θ), 0.0, 0.0, rng)
        fire!(tw, :a, 0.7)
        return pathloglikelihood(tw, 0.7)
    end

    θ0 = 1.3
    @test isapprox(pl_component(θ0, 1), tw_value(θ0); rtol=1e-12, atol=1e-12)
    @test isapprox(ForwardDiff.derivative(θ -> pl_component(θ, 1), θ0),
                   ForwardDiff.derivative(tw_value, θ0); rtol=1e-10, atol=1e-12)

    # The back-compat constructors must still yield Float64 accumulators.
    @test TrajectoryWatcher{Symbol,Float64}() isa TrajectoryWatcher{Symbol,Float64,Float64}
    @test PathLikelihoods{Symbol,Float64}(2) isa PathLikelihoods{Symbol,Float64,Float64}
end
