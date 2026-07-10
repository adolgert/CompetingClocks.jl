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
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t)
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
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    # Path B: sum steploglikelihood over a bare TrackWatcher on the same trace.
    function loglik_step(θ)
        la, lb = θ[1], θ[2]
        track = TrackWatcher{Symbol,Float64}()
        enable!(track, :a, Exponential(1 / la), 0.0, 0.0)
        enable!(track, :b, Exponential(1 / lb), 0.0, 0.0)
        total = zero(eltype(θ))
        t0 = 0.0
        for (t, k) in trace
            total += steploglikelihood(track, t0, t, k)
            fire!(track, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(track, k, d, t, t)
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
        enable!(tw, :w, Weibull(shape, scale), 0.0, 0.0)
        enable!(tw, :e, Exponential(1 / rate), 0.0, 0.0)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :w ? Weibull(shape, scale) : Exponential(1 / rate)
            enable!(tw, k, d, t, t)
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
        enable!(pl, :a, dists, 0.0, 0.0)
        fire!(pl, :a, 0.7)
        return pathloglikelihood(pl, 0.7)[idx]
    end

    function tw_value(θ)
        tw = TrajectoryWatcher{Symbol,Float64,typeof(θ)}()
        enable!(tw, :a, Exponential(1 / θ), 0.0, 0.0)
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


@safetestset autodiff_context_gradient_equals_watcher_and_analytic = "SamplingContext path-likelihood gradient equals the watcher and analytic score" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: TrajectoryWatcher
    using ForwardDiff
    using Random

    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    # The user-facing path: the differentiable distributions reach the likelihood
    # watcher while the sampler receives a primal (value-only) shadow. Using
    # FirstToFireMethod matters: it draws a firing time into a Float64 heap at
    # enable!, so the closure only runs at all if the primal boundary handed the
    # sampler a plain-Float64 distribution.
    function loglik_ctx(θ)
        la, lb = θ[1], θ[2]
        ctx = SamplingContext(SamplerBuilder(Symbol, Float64; path_likelihood=true,
            likelihood_eltype=eltype(θ), method=FirstToFireMethod()), Xoshiro(1))
        enable!(ctx, :a, Exponential(1 / la))
        enable!(ctx, :b, Exponential(1 / lb))
        for (t, k) in trace
            fire!(ctx, k, t)               # advances ctx.time to t
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(ctx, k, d)             # te defaults to now == t, matching the watcher
        end
        return pathloglikelihood(ctx, trace[end][1])
    end

    # The watcher-only path from the existing tests, for a direct comparison.
    rng = Xoshiro(1)
    function loglik_watcher(θ)
        la, lb = θ[1], θ[2]
        tw = TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    θ = [0.7, 0.4]
    g_ctx = ForwardDiff.gradient(loglik_ctx, θ)
    analytic = [3 / 0.7 - 4.5, 2 / 0.4 - 4.5]
    @test isapprox(g_ctx, analytic; rtol=1e-10, atol=1e-12)
    @test isapprox(g_ctx, ForwardDiff.gradient(loglik_watcher, θ); rtol=1e-10, atol=1e-12)
end


@safetestset autodiff_context_weibull_matches_finite_difference = "SamplingContext Weibull-race gradient matches finite differences" begin
    using Distributions
    using CompetingClocks
    using ForwardDiff
    using Random

    # A non-exponential clock exercises logpdf/logccdf derivatives through the
    # context's primal boundary, where the sampler only ever sees Weibull{Float64}.
    trace = [(0.5, :e), (1.2, :w), (2.0, :e)]

    function loglik(θ)
        shape, scale, rate = θ[1], θ[2], θ[3]
        ctx = SamplingContext(SamplerBuilder(Symbol, Float64; path_likelihood=true,
            likelihood_eltype=eltype(θ), method=FirstToFireMethod()), Xoshiro(1))
        enable!(ctx, :w, Weibull(shape, scale))
        enable!(ctx, :e, Exponential(1 / rate))
        for (t, k) in trace
            fire!(ctx, k, t)
            d = k == :w ? Weibull(shape, scale) : Exponential(1 / rate)
            enable!(ctx, k, d)
        end
        return pathloglikelihood(ctx, trace[end][1])
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


@safetestset autodiff_context_hessian_equals_watcher_hessian = "SamplingContext path-likelihood Hessian equals the watcher Hessian" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: TrajectoryWatcher
    using ForwardDiff
    using Random

    # Hessian differentiates a gradient, so the parameters arrive as Dual{Dual}.
    # This proves primal_distribution strips NESTED duals down to Float64 before
    # the sampler sees them.
    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    function loglik_ctx(θ)
        la, lb = θ[1], θ[2]
        ctx = SamplingContext(SamplerBuilder(Symbol, Float64; path_likelihood=true,
            likelihood_eltype=eltype(θ), method=FirstToFireMethod()), Xoshiro(1))
        enable!(ctx, :a, Exponential(1 / la))
        enable!(ctx, :b, Exponential(1 / lb))
        for (t, k) in trace
            fire!(ctx, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(ctx, k, d)
        end
        return pathloglikelihood(ctx, trace[end][1])
    end

    rng = Xoshiro(1)
    function loglik_watcher(θ)
        la, lb = θ[1], θ[2]
        tw = TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
        enable!(tw, :a, Exponential(1 / la), 0.0, 0.0)
        enable!(tw, :b, Exponential(1 / lb), 0.0, 0.0)
        for (t, k) in trace
            fire!(tw, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(tw, k, d, t, t)
        end
        return pathloglikelihood(tw, trace[end][1])
    end

    θ = [0.7, 0.4]
    H_ctx = ForwardDiff.hessian(loglik_ctx, θ)
    @test all(isfinite, H_ctx)
    @test isapprox(H_ctx, ForwardDiff.hessian(loglik_watcher, θ); rtol=1e-10, atol=1e-12)
end


@safetestset autodiff_context_walkup_then_sample_is_finite = "Sampling a continuation at the primal point composes with a dual score" begin
    using Distributions
    using CompetingClocks
    using ForwardDiff
    using Random

    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    # Replay a prefix of the trace, then let the sampler DRAW the continuation.
    # The sampler draws at the primal parameter point (a Float64 time), and that
    # constant time then flows into a Dual-valued path score. The invariant is
    # only that a real, finite gradient survives the round trip.
    function loglik(θ)
        la, lb = θ[1], θ[2]
        ctx = SamplingContext(SamplerBuilder(Symbol, Float64; path_likelihood=true,
            likelihood_eltype=eltype(θ), method=FirstToFireMethod()), Xoshiro(20250706))
        enable!(ctx, :a, Exponential(1 / la))
        enable!(ctx, :b, Exponential(1 / lb))
        for (t, k) in trace[1:3]
            fire!(ctx, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(ctx, k, d)
        end
        when, which = next(ctx)
        fire!(ctx, which, when)
        return pathloglikelihood(ctx, when)
    end

    θ = [0.7, 0.4]
    g = ForwardDiff.gradient(loglik, θ)
    @test all(isfinite, g)
    @test any(!iszero, g)   # the drawn continuation still carries score information
end


@safetestset autodiff_step_likelihood_guard_forces_watcher = "step likelihood keeps derivatives even with a natively-capable sampler" begin
    using Distributions
    using CompetingClocks
    using ForwardDiff
    using Random

    trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]

    # NextReactionMethod builds a sampler that CAN compute steploglikelihood
    # itself. Without the guard, the context would use that native path, which
    # runs on the primal (value-only) distributions and returns an underived
    # Float64 -> a zero gradient. The guard forces a TrackWatcher instead, which
    # holds the Dual-parameterized distributions, so the gradient survives.
    function loglik_step(θ)
        la, lb = θ[1], θ[2]
        ctx = SamplingContext(SamplerBuilder(Symbol, Float64; step_likelihood=true,
            likelihood_eltype=eltype(θ), method=NextReactionMethod()), Xoshiro(1))
        enable!(ctx, :a, Exponential(1 / la))
        enable!(ctx, :b, Exponential(1 / lb))
        total = zero(eltype(θ))
        for (t, k) in trace
            total += steploglikelihood(ctx, t, k)   # ctx.time is the previous fire time
            fire!(ctx, k, t)
            d = k == :a ? Exponential(1 / la) : Exponential(1 / lb)
            enable!(ctx, k, d)
        end
        return total
    end

    θ = [0.7, 0.4]
    g = ForwardDiff.gradient(loglik_step, θ)
    # Summed step likelihoods over the whole trace equal the path likelihood,
    # whose exponential-race score is n_k/λ_k - t_N.
    analytic = [3 / 0.7 - 4.5, 2 / 0.4 - 4.5]
    @test any(!iszero, g)   # a missing guard would zero this out
    @test isapprox(g, analytic; rtol=1e-10, atol=1e-12)
end


@safetestset primal_distribution_strips_dual_parameters = "primal_distribution is identity on numbers and strips AD tracer parameters" begin
    using Distributions
    using CompetingClocks
    using CompetingClocks: primal_distribution
    using ForwardDiff

    # For an ordinary number-typed distribution nothing changes: same object.
    w = Weibull(1.5, 2.0)
    @test primal_distribution(w) === w

    # Dual parameters are stripped to their value parts, giving a Float64 family.
    D = ForwardDiff.Dual
    pw = primal_distribution(Weibull(D(1.5), D(2.0)))
    @test pw isa Weibull{Float64}
    @test params(pw) == (1.5, 2.0)

    # Truncated is a wrapper type: the underlying distribution and any dual bounds
    # are stripped, and open bounds pass through as `nothing`.
    pt = primal_distribution(truncated(Exponential(D(2.0)); lower=D(0.5), upper=3.0))
    @test pt isa Truncated
    @test pt.untruncated isa Exponential{Float64}
    @test pt.lower == 0.5 && pt.upper == 3.0

    # Nested duals (as produced by ForwardDiff.hessian) collapse recursively.
    d1 = ForwardDiff.Dual{:t1}(1.5, 1.0)
    d2 = ForwardDiff.Dual{:t2}(d1, d1)
    pn = primal_distribution(Weibull(d2, d2))
    @test pn isa Weibull{Float64}
end


@safetestset builder_rejects_likelihood_eltype_without_watcher = "SamplerBuilder rejects a differentiable eltype when no likelihood is requested" begin
    using CompetingClocks
    using ForwardDiff

    DualT = typeof(ForwardDiff.Dual(1.0, 1.0))

    # A non-Float64 accumulator is meaningless without a likelihood watcher.
    @test_throws ArgumentError SamplerBuilder(Symbol, Float64; likelihood_eltype=DualT)

    # It is accepted with any of the likelihood options.
    @test SamplerBuilder(Symbol, Float64; path_likelihood=true, likelihood_eltype=DualT) isa SamplerBuilder
    @test SamplerBuilder(Symbol, Float64; step_likelihood=true, likelihood_eltype=DualT) isa SamplerBuilder
    @test SamplerBuilder(Symbol, Float64; likelihood_cnt=2, likelihood_eltype=DualT) isa SamplerBuilder
end
