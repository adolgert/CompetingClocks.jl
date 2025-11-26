# HMC over Event Paths with Gen.jl

This guide shows how to perform Hamiltonian Monte Carlo (HMC) inference over event times in a continuous-time discrete-event system (CTDES), integrating CompetingClocks with [Gen.jl](https://www.gen.dev/) and [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).

## Statistical Setup

Consider a CTDES with clock keys ``k \in \mathcal{K}`` where each enabled clock has a waiting-time distribution ``F_k(\cdot; \theta)``. A trajectory is a sequence of events:

```math
\text{path} = \big((k_1, t_1), (k_2, t_2), \ldots, (k_M, t_M)\big), \quad 0 < t_1 < \cdots < t_M \le T
```

Given observed data ``y`` derived from the path (e.g., a noisy completion time), the posterior over event times and parameters is:

```math
\log p(\theta, \text{path} \mid y) = \log p_\theta(\text{path}) + \log p(y \mid \text{path}, \theta) + \log p(\theta) + C
```

HMC treats the event times as a continuous vector ``q``, introduces momenta ``p``, and simulates Hamiltonian dynamics using gradients ``\nabla_q \log p(\theta, \text{path} \mid y)``.

Two integration strategies with Gen.jl:

1. **External HMC over event times**: Gen handles discrete structure and priors; a separate HMC kernel updates event times.
2. **Gen-native gradients**: Wrap the path log-likelihood as a custom Gen distribution with `logpdf_grad`.

This guide demonstrates approach (1).

## Example Model: Two-Step Reaction A → B → C

A minimal chain reaction:
- At ``t=0``, the system is in state `:A`
- Event `:ab` fires once with ``\text{Weibull}(k_{ab}, \lambda_{ab})`` waiting time
- After `:ab` fires, event `:bc` is enabled with ``\text{Weibull}(k_{bc}, \lambda_{bc})`` waiting time
- We observe a noisy completion time ``y \sim \mathcal{N}(t_2, \sigma_{\text{obs}}^2)``

The goal is posterior inference over ``(t_1, t_2)`` given ``y``.

## Dependencies

```julia
using Random
using Distributions
using CompetingClocks
using ForwardDiff
using LinearAlgebra
using Gen
```

## Model Definitions

```julia
# Clock key for the two-step reaction
struct ABClock
    kind::Symbol  # :ab or :bc
end

# Weibull parameters for both transitions
struct ABParams{T}
    k_ab::T
    lambda_ab::T
    k_bc::T
    lambda_bc::T
end

# Model state
mutable struct ABModel{T}
    t::T
    state::Symbol
    params::ABParams{T}
end
```

## Event Handling

```julia
function initialize_ab!(model::ABModel{T}, sampler) where {T}
    model.t = zero(T)
    model.state = :A
    p = model.params
    enable!(sampler, ABClock(:ab), Weibull(p.k_ab, p.lambda_ab))
    return nothing
end

function handle_ab_event!(model::ABModel{T}, sampler, which::ABClock, when::T) where {T}
    model.t = when
    p = model.params

    if which.kind == :ab
        model.state = :B
        enable!(sampler, ABClock(:bc), Weibull(p.k_bc, p.lambda_bc))
    elseif which.kind == :bc
        model.state = :C
    else
        error("Unknown clock kind $(which.kind)")
    end
end
```

## Path Log-Likelihood

For this simple model, the path log-likelihood has a closed form:

```math
\log p(t_1, t_2 \mid \theta) = \log f_{ab}(t_1) + \log f_{bc}(t_2 - t_1)
```

where ``f_{ab}`` and ``f_{bc}`` are the Weibull densities.

```julia
function log_path_ab(events::AbstractVector, params::ABParams)
    @assert length(events) == 2
    t1 = events[1].time
    t2 = events[2].time

    dist_ab = Weibull(params.k_ab, params.lambda_ab)
    dist_bc = Weibull(params.k_bc, params.lambda_bc)

    # Use Distributions.logpdf to avoid ambiguity with Gen.logpdf
    return Distributions.logpdf(dist_ab, t1) + Distributions.logpdf(dist_bc, t2 - t1)
end
```

For validation, CompetingClocks can compute the same quantity via replay:

```julia
function log_path_ab_cclock(events::AbstractVector, params::ABParams{Float64})
    rng = Xoshiro(0)
    sampler = SamplingContext(ABClock, Float64, rng; path_likelihood=true)

    model = ABModel(0.0, :A, params)
    initialize_ab!(model, sampler)

    for e in events
        fire!(sampler, e.evt, Float64(e.time))
        handle_ab_event!(model, sampler, e.evt, Float64(e.time))
    end

    return pathloglikelihood(sampler, Float64(events[end].time))
end
```

!!! note "AD Compatibility"
    The `TrajectoryWatcher` inside `SamplingContext` currently stores `loglikelihood::Float64`, which prevents direct use with ForwardDiff dual numbers. For HMC, compute the path log-likelihood manually as shown in `log_path_ab`.

## Unconstrained Parameterization

To maintain ordered, positive event times during HMC, map unconstrained ``u \in \mathbb{R}^2`` to times:

```math
t_1 = e^{u_1}, \quad t_2 = t_1 + e^{u_2}
```

with log-Jacobian ``\log |J| = u_1 + u_2``.

```julia
function times_from_u(u::AbstractVector)
    @assert length(u) == 2
    t1 = exp(u[1])
    dt = exp(u[2])
    t2 = t1 + dt
    logjac = u[1] + u[2]
    return (t1, t2, logjac)
end
```

## Log Posterior

```julia
function logpost_u(u::AbstractVector, params::ABParams, y_obs::Float64, sigma_obs::Float64)
    t1, t2, logjac = times_from_u(u)

    events = [
        (evt=ABClock(:ab), time=t1),
        (evt=ABClock(:bc), time=t2)
    ]

    lp_path = log_path_ab(events, params)
    lp_obs = Distributions.logpdf(Normal(t2, sigma_obs), y_obs)

    return lp_path + lp_obs + logjac
end

# Callable wrapper for ForwardDiff
struct PosteriorWrapper
    params::ABParams{Float64}
    y_obs::Float64
    sigma_obs::Float64
end

(pw::PosteriorWrapper)(u) = logpost_u(u, pw.params, pw.y_obs, pw.sigma_obs)
```

## HMC Implementation

A basic leapfrog integrator with Metropolis accept/reject:

```julia
function hmc_step(u_current::Vector{Float64}, logpost; epsilon::Float64=0.02, L::Int=10)
    q = copy(u_current)
    p = randn(length(q))
    current_p = copy(p)

    grad_logpost(x) = ForwardDiff.gradient(logpost, x)

    # Leapfrog integration
    grad = grad_logpost(q)
    p .+= (epsilon / 2.0) .* grad

    for l in 1:L
        q .+= epsilon .* p
        grad = grad_logpost(q)
        if l != L
            p .+= epsilon .* grad
        end
    end

    p .+= (epsilon / 2.0) .* grad
    p .= -p

    # Metropolis accept/reject
    H(q, p) = -logpost(q) + 0.5 * dot(p, p)
    current_H = H(u_current, current_p)
    proposed_H = H(q, p)

    if log(rand()) < current_H - proposed_H
        return q, true
    else
        return u_current, false
    end
end
```

## Running HMC

```julia
function demo_hmc()
    true_params = ABParams(1.5, 3.0, 2.0, 1.5)

    # Generate synthetic data
    rng = MersenneTwister(42)
    t1_true = rand(rng, Weibull(true_params.k_ab, true_params.lambda_ab))
    t2_true = t1_true + rand(rng, Weibull(true_params.k_bc, true_params.lambda_bc))
    sigma_obs = 0.1
    y_obs = rand(rng, Normal(t2_true, sigma_obs))

    pw = PosteriorWrapper(true_params, y_obs, sigma_obs)
    u = [log(0.5 * y_obs), log(0.5 * y_obs)]

    samples = Vector{Vector{Float64}}()
    acc = 0
    n_iter = 500

    for iter in 1:n_iter
        u, accepted = hmc_step(u, pw; epsilon=0.01, L=20)
        acc += accepted ? 1 : 0
        push!(samples, copy(u))
    end

    # Posterior summary
    burn_in = 100
    t1_samples = [times_from_u(s)[1] for s in samples[burn_in+1:end]]
    t2_samples = [times_from_u(s)[2] for s in samples[burn_in+1:end]]

    println("Acceptance rate: $(acc / n_iter)")
    println("Posterior mean t1: $(mean(t1_samples)) (true: $t1_true)")
    println("Posterior mean t2: $(mean(t2_samples)) (true: $t2_true)")

    return samples
end
```

## Gen.jl Integration

To infer model parameters ``\theta`` while marginalizing over event times, use HMC as an inner loop within Gen's inference framework.

### Factor Distribution

Gen.jl doesn't export a `factor` function for adding log-weights. Define a custom distribution:

```julia
struct FactorDist <: Gen.Distribution{Nothing} end
const factor_dist = FactorDist()

Gen.random(::FactorDist, logweight::Float64) = nothing
Gen.logpdf(::FactorDist, ::Nothing, logweight::Float64) = logweight
Gen.is_discrete(::FactorDist) = true
Gen.has_output_grad(::FactorDist) = false
Gen.has_argument_grads(::FactorDist) = (false,)
```

### Approximate Marginal Likelihood

Use a short HMC run to approximate ``\log p(y \mid \theta)``:

```julia
function loglik_given_params(params::ABParams, y_obs::Float64, sigma_obs::Float64)
    pw = PosteriorWrapper(params, y_obs, sigma_obs)
    u = [log(0.5 * y_obs), log(0.5 * y_obs)]

    for _ in 1:50
        u, _ = hmc_step(u, pw; epsilon=0.01, L=10)
    end

    return pw(u)
end
```

### Gen Model

```julia
@gen function ab_params_model()
    delta_k_ab = @trace(normal(0.0, 0.1), :delta_k_ab)
    delta_lambda_ab = @trace(normal(0.0, 0.1), :delta_lambda_ab)
    delta_k_bc = @trace(normal(0.0, 0.1), :delta_k_bc)
    delta_lambda_bc = @trace(normal(0.0, 0.1), :delta_lambda_bc)

    base = ABParams(1.5, 3.0, 2.0, 1.5)
    return ABParams(
        base.k_ab + delta_k_ab,
        base.lambda_ab + delta_lambda_ab,
        base.k_bc + delta_k_bc,
        base.lambda_bc + delta_lambda_bc
    )
end

@gen function joint_model(y_obs::Float64, sigma_obs::Float64)
    params = @trace(ab_params_model(), :params)
    lp_y = loglik_given_params(params, y_obs, sigma_obs)
    @trace(factor_dist(lp_y), :factor)
    return params
end
```

### MCMC over Parameters

```julia
function run_gen_with_hmc(y_obs, sigma_obs)
    tr, _ = generate(joint_model, (y_obs, sigma_obs))

    @gen function proposal_kernel(tr)
        {:params => :delta_k_ab} ~ normal(tr[:params => :delta_k_ab], 0.05)
        {:params => :delta_lambda_ab} ~ normal(tr[:params => :delta_lambda_ab], 0.05)
        {:params => :delta_k_bc} ~ normal(tr[:params => :delta_k_bc], 0.05)
        {:params => :delta_lambda_bc} ~ normal(tr[:params => :delta_lambda_bc], 0.05)
        return nothing
    end

    n_steps = 100
    for i in 1:n_steps
        tr, _ = Gen.mh(tr, proposal_kernel, ())
    end

    return tr
end
```

## Summary

This approach combines:

- **CompetingClocks** for defining the CTDES and computing path log-likelihoods
- **ForwardDiff** for gradients through the log-posterior
- **HMC** for efficient sampling of continuous event times
- **Gen.jl** for probabilistic programming over parameters

The key constraint is that path log-likelihood computation must be AD-compatible. For simple models, compute it analytically. For complex models, consider modifying CompetingClocks' `TrajectoryWatcher` to support generic numeric types.

## Full Example

See [`examples/gen_hmc.jl`](https://github.com/adolgert/CompetingClocks.jl/blob/main/examples/gen_hmc.jl) for a complete, runnable version.
