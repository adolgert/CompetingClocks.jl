# Gen.jl + CompetingClocks + HMC for Bayesian inference over event paths
#
# This example demonstrates Hamiltonian Monte Carlo (HMC) over event times
# in a continuous-time discrete-event system (CTDES), combining:
#   - CompetingClocks for path log-likelihood computation
#   - ForwardDiff for automatic differentiation of the log-posterior
#   - Gen.jl for parameter inference
#
# The model is a simple two-step reaction: A -> B -> C
#   - Event :ab fires once (A->B), Weibull waiting time
#   - Event :bc fires once (B->C), Weibull waiting time
#   - We observe a noisy completion time y ~ Normal(t2, sigma_obs)
#
# Run with: julia --project=examples examples/gen_hmc.jl

println("Loading packages...")
using Random
using Distributions
using CompetingClocks
using ForwardDiff
using LinearAlgebra
println("Packages loaded.")

# ===========================================================================
# Part 1: Clock key and model state definitions
# ===========================================================================

println("=" ^ 60)
println("Part 1: Model definitions")
println("=" ^ 60)

"""
    ABClock

Clock key for the two-step reaction A -> B -> C.
"""
struct ABClock
    kind::Symbol  # :ab or :bc
end

"""
    ABParams{T}

Parameters for the two Weibull waiting times.
- A -> B ~ Weibull(k_ab, lambda_ab)
- B -> C ~ Weibull(k_bc, lambda_bc)
"""
struct ABParams{T}
    k_ab::T
    lambda_ab::T
    k_bc::T
    lambda_bc::T
end

"""
    ABModel{T}

Minimal state for the A->B->C system.
"""
mutable struct ABModel{T}
    t::T
    state::Symbol
    params::ABParams{T}
end

println("Defined ABClock, ABParams, ABModel types.")

# ===========================================================================
# Part 2: Enabling and handling events
# ===========================================================================

println("=" ^ 60)
println("Part 2: Event handling functions")
println("=" ^ 60)

"""
    initialize_ab!(model, sampler)

At t = 0, system is in state :A and only the :ab clock is enabled.
"""
function initialize_ab!(model::ABModel{T}, sampler) where {T}
    model.t = zero(T)
    model.state = :A
    p = model.params
    dist_ab = Weibull(p.k_ab, p.lambda_ab)
    # Enable A->B starting now
    enable!(sampler, ABClock(:ab), dist_ab)
    return nothing
end

"""
    handle_ab_event!(model, sampler, which, when)

Update state and enabled clocks for the two events.
Assumes `fire!(sampler, which, when)` has just been called.
"""
function handle_ab_event!(
    model::ABModel{T},
    sampler,
    which::ABClock,
    when::T
) where {T}

    model.t = when
    p = model.params

    if which.kind == :ab
        # Transition A->B: disable :ab, enable :bc
        model.state = :B
        # `fire!` has already disabled :ab in the sampler.
        dist_bc = Weibull(p.k_bc, p.lambda_bc)
        enable!(sampler, ABClock(:bc), dist_bc)

    elseif which.kind == :bc
        # Transition B->C: final state, no further clocks
        model.state = :C
        # `fire!` already disabled :bc.

    else
        error("Unknown clock kind $(which.kind)")
    end
end

println("Defined initialize_ab! and handle_ab_event!.")

# ===========================================================================
# Part 3: Path log-likelihood from an event list
# ===========================================================================

println("=" ^ 60)
println("Part 3: Path log-likelihood computation")
println("=" ^ 60)

"""
    log_path_ab(events, params)

Compute path log-likelihood log p(path | params) for the given ordered
event list `events`, where each element is a NamedTuple: (evt=ABClock, time=T).

This is computed manually to be AD-friendly with ForwardDiff.
For the A->B->C model:
  - Event :ab fires at t1, Weibull(k_ab, lambda_ab) waiting time from t=0
  - Event :bc fires at t2, Weibull(k_bc, lambda_bc) waiting time from t=t1

The log-likelihood is:
  log p(t1, t2 | params) = logpdf(Weibull(k_ab, lambda_ab), t1)
                         + logpdf(Weibull(k_bc, lambda_bc), t2 - t1)
"""
function log_path_ab(events::AbstractVector, params::ABParams)
    @assert length(events) == 2

    t1 = events[1].time  # Time of :ab event
    t2 = events[2].time  # Time of :bc event

    # Weibull distributions for the two events
    dist_ab = Weibull(params.k_ab, params.lambda_ab)
    dist_bc = Weibull(params.k_bc, params.lambda_bc)

    # Log-likelihood: logpdf of each waiting time
    # Use Distributions.logpdf to avoid ambiguity with Gen.logpdf
    ll = Distributions.logpdf(dist_ab, t1) + Distributions.logpdf(dist_bc, t2 - t1)

    return ll
end

"""
    log_path_ab_cclock(events, params)

Alternative version using CompetingClocks' TrajectoryWatcher directly.
Note: This does NOT work with ForwardDiff because TrajectoryWatcher
has Float64 fields internally. Use for validation with Float64 only.
"""
function log_path_ab_cclock(events::AbstractVector, params::ABParams{Float64})
    @assert !isempty(events)

    # Deterministic RNG seed so that log-density is pure
    rng = Xoshiro(0)

    # path_likelihood=true tells CompetingClocks to track path log-likelihood
    sampler = SamplingContext(ABClock, Float64, rng; path_likelihood=true)

    # Initialize model and first clock
    model = ABModel(
        0.0,
        :A,
        params
    )
    initialize_ab!(model, sampler)

    # Replay events in order
    for e in events
        which = e.evt
        when = Float64(e.time)  # Convert to Float64
        # Advance sampler time and notify it of which clock fired
        fire!(sampler, which, when)
        # Update model and (de)enable clocks
        handle_ab_event!(model, sampler, which, when)
    end

    t_end = Float64(events[end].time)
    return pathloglikelihood(sampler, t_end)
end

# Quick test of log_path_ab
test_params = ABParams(1.5, 3.0, 2.0, 1.5)
test_events = [
    (evt=ABClock(:ab), time=1.0),
    (evt=ABClock(:bc), time=2.5)
]
test_logp_manual = log_path_ab(test_events, test_params)
test_logp_cclock = log_path_ab_cclock(test_events, test_params)
println("Manual log-likelihood:       $test_logp_manual")
println("CompetingClocks log-likelihood: $test_logp_cclock")
println("Match: $(isapprox(test_logp_manual, test_logp_cclock))")
println()

# ===========================================================================
# Part 4: Log posterior over event times and HMC
# ===========================================================================

println("=" ^ 60)
println("Part 4: Log posterior and HMC implementation")
println("=" ^ 60)

"""
    times_from_u(u)

Map unconstrained vector u = [u1, u2] to ordered times t1, t2:
    t1 = exp(u1)
    t2 = t1 + exp(u2)

Returns (t1, t2, logjac) where logjac is the log-Jacobian of this transform.
"""
function times_from_u(u::AbstractVector)
    @assert length(u) == 2
    t1 = exp(u[1])
    dt = exp(u[2])
    t2 = t1 + dt
    logjac = u[1] + u[2]  # log|det J|
    return (t1, t2, logjac)
end

"""
    logpost_u(u, params, y_obs, sigma_obs)

Log posterior (up to a constant) as a function of unconstrained u:
    log p(u | y) = log p(path(u)) + log p(y | t2) + logJacobian + const.

Here:
- path(u) is the two-event path with times t1, t2.
- prior over u is taken to be flat in this example.
"""
function logpost_u(
    u::AbstractVector,
    params::ABParams,
    y_obs::Float64,
    sigma_obs::Float64
)
    t1, t2, logjac = times_from_u(u)

    # Build event list expected by log_path_ab
    events = [
        (evt=ABClock(:ab), time=t1),
        (evt=ABClock(:bc), time=t2)
    ]

    # Path log-likelihood under the CTDES
    lp_path = log_path_ab(events, params)

    # Observation: noisy measurement of completion time t2
    lp_obs = Distributions.logpdf(Normal(t2, sigma_obs), y_obs)

    # Flat prior over u: lp_prior = 0.0
    return lp_path + lp_obs + logjac
end

"""
    PosteriorWrapper

Convenience wrapper for ForwardDiff gradient computation.
"""
struct PosteriorWrapper
    params::ABParams{Float64}
    y_obs::Float64
    sigma_obs::Float64
end

(pw::PosteriorWrapper)(u) = logpost_u(u, pw.params, pw.y_obs, pw.sigma_obs)

"""
    hmc_step(u_current, logpost; epsilon, L)

One HMC step on u using leapfrog integrator.
Returns (u_new, accepted::Bool).
"""
function hmc_step(
    u_current::Vector{Float64},
    logpost;  # Any callable (Function or functor)
    epsilon::Float64=0.02,
    L::Int=10
)
    # Current position and momentum
    q = copy(u_current)
    p = randn(length(q))
    current_p = copy(p)

    # Helper: gradient of logpost
    grad_logpost(x) = ForwardDiff.gradient(logpost, x)

    # Initial gradient
    grad = grad_logpost(q)

    # Half step in momentum
    p .+= (epsilon / 2.0) .* grad

    # Full leapfrog steps
    for l in 1:L
        # Full step in position
        q .+= epsilon .* p
        # Gradient at new position
        grad = grad_logpost(q)
        if l != L
            p .+= epsilon .* grad
        end
    end

    # Final half step in momentum
    p .+= (epsilon / 2.0) .* grad

    # Negate momentum for symmetry
    p .= -p

    # Hamiltonians
    function H(q, p)
        U = -logpost(q)
        K = 0.5 * dot(p, p)
        return U + K
    end

    current_H = H(u_current, current_p)
    proposed_H = H(q, p)

    accept_logprob = current_H - proposed_H
    if log(rand()) < accept_logprob
        return q, true
    else
        return u_current, false
    end
end

println("Defined times_from_u, logpost_u, PosteriorWrapper, hmc_step.")
println()

# ===========================================================================
# Part 5: Running a small HMC chain
# ===========================================================================

println("=" ^ 60)
println("Part 5: Demo HMC chain")
println("=" ^ 60)

function demo_hmc()
    # "True" parameters for data generation
    true_params = ABParams(1.5, 3.0, 2.0, 1.5)

    # Generate one synthetic path and noisy completion time
    rng = MersenneTwister(42)
    t1_true = rand(rng, Weibull(true_params.k_ab, true_params.lambda_ab))
    t2_true = t1_true + rand(rng, Weibull(true_params.k_bc, true_params.lambda_bc))
    sigma_obs = 0.1
    y_obs = rand(rng, Normal(t2_true, sigma_obs))

    println("True times: t1 = $(round(t1_true, digits=4)), t2 = $(round(t2_true, digits=4))")
    println("Observed y = $(round(y_obs, digits=4))")

    # Use same params for inference
    infer_params = ABParams(1.5, 3.0, 2.0, 1.5)

    pw = PosteriorWrapper(infer_params, y_obs, sigma_obs)

    # Initial u (log-times) guess
    u = [log(0.5 * y_obs), log(0.5 * y_obs)]

    samples = Vector{Vector{Float64}}()
    acc = 0
    n_iter = 500  # Reduced for faster demo

    println("Running HMC chain with $n_iter iterations...")

    for iter in 1:n_iter
        u, accepted = hmc_step(u, pw; epsilon=0.01, L=20)
        acc += accepted ? 1 : 0
        push!(samples, copy(u))
    end

    acceptance_rate = acc / n_iter
    println("Acceptance rate = $(round(acceptance_rate, digits=3))")

    # Map samples to times
    t1_samples = [times_from_u(s)[1] for s in samples]
    t2_samples = [times_from_u(s)[2] for s in samples]

    # Compute posterior means (burn-in first 100)
    burn_in = min(100, n_iter ÷ 5)
    t1_mean = mean(t1_samples[burn_in+1:end])
    t2_mean = mean(t2_samples[burn_in+1:end])

    println("Posterior mean t1 = $(round(t1_mean, digits=4)) (true: $(round(t1_true, digits=4)))")
    println("Posterior mean t2 = $(round(t2_mean, digits=4)) (true: $(round(t2_true, digits=4)))")

    return samples, t1_true, t2_true
end

samples_hmc, t1_true, t2_true = demo_hmc()
println()

# ===========================================================================
# Part 6: Combining with Gen.jl
# ===========================================================================

println("Loading Gen...")
using Gen
println("Gen loaded.")

println("=" ^ 60)
println("Part 6: Gen.jl integration")
println("=" ^ 60)

"""
    loglik_given_params(params, y_obs, sigma_obs)

For fixed params, uses a short HMC run to approximate one sample of the
latent path and returns its log posterior (up to const).

In a real application you might:
- run a longer chain and average,
- or use MAP event times via optimization
to get a better approximation to log p(y | params).
"""
function loglik_given_params(params::ABParams, y_obs::Float64, sigma_obs::Float64)
    pw = PosteriorWrapper(params, y_obs, sigma_obs)
    u = [log(0.5 * y_obs), log(0.5 * y_obs)]

    # A small number of HMC steps to move toward high-density region
    for _ in 1:50
        u, _ = hmc_step(u, pw; epsilon=0.01, L=10)
    end

    # Use final state as Laplace-like approximation to integral over paths
    return pw(u)
end

# ===========================================================================
# Custom "Factor" distribution for adding log-weights in Gen
# ===========================================================================

"""
    FactorDist

A custom Gen distribution that adds a deterministic log-weight (factor) to the trace.
The "value" is always nothing, and logpdf returns the specified weight.
This is a workaround since Gen.jl doesn't have a built-in `factor` function.
"""
struct FactorDist <: Gen.Distribution{Nothing} end
const factor_dist = FactorDist()

function Gen.random(::FactorDist, logweight::Float64)
    return nothing
end

function Gen.logpdf(::FactorDist, value::Nothing, logweight::Float64)
    return logweight
end

Gen.is_discrete(::FactorDist) = true
Gen.has_output_grad(::FactorDist) = false
Gen.has_argument_grads(::FactorDist) = (false,)

# ===========================================================================
# Gen models
# ===========================================================================

# Prior over parameters (sampling small perturbations around nominal values)
@gen function ab_params_model()
    # Sample perturbations around nominal parameter values
    delta_k_ab = @trace(normal(0.0, 0.1), :delta_k_ab)
    delta_lambda_ab = @trace(normal(0.0, 0.1), :delta_lambda_ab)
    delta_k_bc = @trace(normal(0.0, 0.1), :delta_k_bc)
    delta_lambda_bc = @trace(normal(0.0, 0.1), :delta_lambda_bc)

    base = ABParams(1.5, 3.0, 2.0, 1.5)
    params = ABParams(
        base.k_ab + delta_k_ab,
        base.lambda_ab + delta_lambda_ab,
        base.k_bc + delta_k_bc,
        base.lambda_bc + delta_lambda_bc
    )
    return params
end

@gen function joint_model(y_obs::Float64, sigma_obs::Float64)
    # Sample parameters under the prior
    params = @trace(ab_params_model(), :params)

    # Approximate marginal log-likelihood of y_obs under those params
    lp_y = loglik_given_params(params, y_obs, sigma_obs)

    # Incorporate this as a factor using our custom distribution
    # This adds lp_y to the trace's log-probability without a random choice
    @trace(factor_dist(lp_y), :factor)

    # Return parameters
    return params
end

# ===========================================================================
# Test Part 6: Running Gen inference over parameters
# ===========================================================================

println("Running Gen inference with HMC-based likelihood...")

function run_gen_with_hmc()
    # Generate synthetic data
    true_params = ABParams(1.5, 3.0, 2.0, 1.5)
    rng = MersenneTwister(123)
    t1_true = rand(rng, Weibull(true_params.k_ab, true_params.lambda_ab))
    t2_true = t1_true + rand(rng, Weibull(true_params.k_bc, true_params.lambda_bc))
    sigma_obs = 0.1
    y_obs = rand(rng, Normal(t2_true, sigma_obs))

    println("Data for Gen model: t2_true = $(round(t2_true, digits=4)), y_obs = $(round(y_obs, digits=4))")

    # Initial trace
    tr, _ = generate(joint_model, (y_obs, sigma_obs))

    # Simple single-site MH kernel on delta parameters
    @gen function proposal_kernel(tr)
        # Propose new deltas with small random walk
        delta_k_ab_p = {:params => :delta_k_ab} ~ normal(tr[:params => :delta_k_ab], 0.05)
        delta_lambda_ab_p = {:params => :delta_lambda_ab} ~ normal(tr[:params => :delta_lambda_ab], 0.05)
        delta_k_bc_p = {:params => :delta_k_bc} ~ normal(tr[:params => :delta_k_bc], 0.05)
        delta_lambda_bc_p = {:params => :delta_lambda_bc} ~ normal(tr[:params => :delta_lambda_bc], 0.05)
        return nothing
    end

    # Basic MH loop
    n_steps = 50  # Reduced for faster demo
    traces = Vector{Any}(undef, n_steps)
    traces[1] = tr

    println("Running MH with $n_steps steps...")
    for i in 2:n_steps
        tr, _ = Gen.mh(tr, proposal_kernel, ())
        traces[i] = tr
        if i % 10 == 0
            params = get_retval(tr)
            println("  Step $i: k_ab=$(round(params.k_ab, digits=3)), lambda_ab=$(round(params.lambda_ab, digits=3))")
        end
    end

    return traces
end

traces = run_gen_with_hmc()

println()
println("=" ^ 60)
println("Example complete!")
println("=" ^ 60)
