# Observation Likelihood for Event Data

This page describes how to compute the likelihood of observed event sequences
using `pathloglikelihood`, and how to integrate this with probabilistic
programming frameworks like [Gen.jl](https://www.gen.dev/) for Bayesian inference.

## Statistical Background

Consider a continuous-time system where multiple events (clocks) compete to fire:

- At time ``t``, a set of clocks ``\mathcal{E}_t`` is enabled.
- Each enabled clock ``e`` has hazard ``\lambda_e(t \mid \mathcal{H}_t, \theta)``,
  where ``\mathcal{H}_t`` is the history up to time ``t`` and ``\theta`` are parameters.
- When a clock fires, the system state changes and enabled clocks are updated.

Given an observed event path over ``[0, T]``,

```math
\mathcal{D} = \{(t_1, e_1), \ldots, (t_n, e_n)\}, \quad 0 < t_1 < \cdots < t_n \le T,
```

the log-likelihood has the standard competing risks form:

```math
\log p(\mathcal{D}\mid\theta) =
\sum_{i=1}^{n} \log \lambda_{e_i}(t_i \mid \mathcal{H}_{t_i^-}, \theta)
- \int_0^T \Lambda(t \mid \mathcal{H}_t, \theta)\,dt
```

where ``\Lambda(t) = \sum_{e \in \mathcal{E}_t} \lambda_e(t)`` is the total hazard.

The first term scores each observed event; the second term accounts for the
probability that no other events occurred. For non-exponential distributions
(Weibull, Gamma, etc.) this integral requires tracking cumulative hazards
from each clock's enabling time.

## CompetingClocks Approach

CompetingClocks maintains these integrals internally. To compute the likelihood:

1. Create a `SamplingContext` with `path_likelihood=true`
2. **Replay** the observed events: call `enable!` and `fire!` in sequence
3. Call `pathloglikelihood(sampler, end_time)` to get ``\log p(\mathcal{D}\mid\theta)``

This avoids manual hazard calculus entirely.

## Example: Single-Machine Reliability

A machine alternates between working (`:up`) and failed (`:down`):

- Time to failure follows ``\text{Weibull}(k, \eta)``
- Time to repair follows ``\text{Exponential}(\mu)``

We observe failure/repair events and want to infer ``\theta = (k, \eta, \mu)``.

### Data Representation

```julia
using Random
using Distributions
using CompetingClocks

const RelEvent = NamedTuple{(:key, :time), Tuple{Symbol, Float64}}

observed_events = RelEvent[
    (key = :fail,   time = 1.2),
    (key = :repair, time = 2.0),
    (key = :fail,   time = 4.5),
    (key = :repair, time = 5.1),
]
obs_end_time = 8.0
```

### Model State

```julia
struct ReliabilityParams
    k_fail::Float64     # Weibull shape
    η_fail::Float64     # Weibull scale
    μ_repair::Float64   # Exponential mean (1/rate)
end

mutable struct ReliabilityModel
    up::Bool
    θ::ReliabilityParams
end

init_model(θ::ReliabilityParams) = ReliabilityModel(true, θ)
```

### Clock Management

```julia
function enable_current_clock!(model::ReliabilityModel, sampler)
    if model.up
        enable!(sampler, :fail, Weibull(model.θ.k_fail, model.θ.η_fail))
    else
        enable!(sampler, :repair, Exponential(model.θ.μ_repair))
    end
end

function step_reliability!(model::ReliabilityModel, sampler, key::Symbol, when::Float64)
    fire!(sampler, key, when)
    model.up = (key == :repair)
    enable_current_clock!(model, sampler)
end
```

### Computing the Log-Likelihood

```julia
function reliability_loglikelihood(θ::ReliabilityParams,
                                   events::Vector{RelEvent},
                                   end_time::Float64)
    sampler = SamplingContext(Symbol, Float64, Xoshiro(1);
                              path_likelihood = true)
    model = init_model(θ)
    enable_current_clock!(model, sampler)

    for evt in events
        step_reliability!(model, sampler, evt.key, evt.time)
    end

    return pathloglikelihood(sampler, end_time)
end
```

Test it:

```julia
θ_true = ReliabilityParams(1.5, 2.0, 1.0)
ll = reliability_loglikelihood(θ_true, observed_events, obs_end_time)
# ll ≈ -5.73
```

## Integration with Gen.jl

Gen.jl provides probabilistic programming with MCMC inference. To use
CompetingClocks likelihoods in Gen, create a custom distribution that
contributes the log-likelihood as a "factor" in the trace.

### Likelihood Factor Distribution

Since Gen's `@factor` macro may not be available in all versions, we define
a custom distribution that adds an arbitrary log-probability to the trace:

```julia
using Gen
using Statistics: mean, std

struct LikelihoodFactor <: Gen.Distribution{Nothing} end
const likelihood_factor = LikelihoodFactor()

Gen.random(::LikelihoodFactor, logpdf_val::Float64) = nothing
Gen.logpdf(::LikelihoodFactor, ::Nothing, logpdf_val::Float64) = logpdf_val
Gen.is_discrete(::LikelihoodFactor) = true
Gen.has_output_grad(::LikelihoodFactor) = false
Gen.has_argument_grads(::LikelihoodFactor) = (false,)
```

### Gen Model

```julia
@gen function reliability_model(events::Vector{RelEvent}, end_time::Float64)
    # Priors on log-parameters (positive support)
    log_k ~ normal(0.0, 0.5)
    log_η ~ normal(0.0, 1.0)
    log_μ ~ normal(0.0, 1.0)

    θ = ReliabilityParams(exp(log_k), exp(log_η), exp(log_μ))

    # Add CompetingClocks likelihood to the trace
    loglike = reliability_loglikelihood(θ, events, end_time)
    {:likelihood} ~ likelihood_factor(loglike)

    return θ
end
```

### Inference

Generate an initial trace:

```julia
(trace, _) = generate(reliability_model, (observed_events, obs_end_time))
println("Initial params: ", get_retval(trace))
println("Log-joint: ", get_score(trace))
```

Run Metropolis-Hastings with Gaussian drift proposals:

```julia
addrs = [:log_k, :log_η, :log_μ]

function run_mh(initial_trace, n_iters; drift_std=0.1)
    samples = Vector{ReliabilityParams}(undef, n_iters)
    current_trace = initial_trace
    accepted = 0

    for i in 1:n_iters
        addr = addrs[(i - 1) % 3 + 1]
        current_val = current_trace[addr]
        proposed_val = current_val + randn() * drift_std
        constraints = choicemap((addr, proposed_val))

        new_trace, weight, _, _ = update(
            current_trace, get_args(current_trace),
            (NoChange(), NoChange()), constraints
        )

        if log(rand()) < weight
            current_trace = new_trace
            accepted += 1
        end
        samples[i] = get_retval(current_trace)
    end

    return samples, accepted / n_iters
end

samples, acc_rate = run_mh(trace, 1000)
println("Acceptance rate: $acc_rate")

# Posterior summary (after burn-in)
post = samples[501:end]
println("k: mean=$(mean(s.k_fail for s in post))")
println("η: mean=$(mean(s.η_fail for s in post))")
println("μ: mean=$(mean(s.μ_repair for s in post))")
```

## Key Points

- `pathloglikelihood` computes the full path likelihood including the
  probability of no events in any remaining observation window.

- The "replay" pattern works for any model structure: simply call `enable!`
  and `fire!` to mirror what would happen during forward simulation.

- The likelihood function is pure Julia and can be used outside Gen for
  maximum likelihood estimation, importance sampling, or other inference methods.

- For models with many parameters or complex posteriors, consider gradient-based
  methods or more sophisticated MCMC schemes.

## See Also

- [Gen.jl Documentation](https://www.gen.dev/docs/)
- [Adding Custom Distributions to Gen](https://www.gen.dev/docs/stable/how_to/custom_distributions/)
- [`examples/gen_observation.jl`](https://github.com/adolgert/CompetingClocks.jl/blob/main/examples/gen_observation.jl) — Complete working example
