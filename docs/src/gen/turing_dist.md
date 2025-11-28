# Bayesian Inference with Turing.jl

This page demonstrates how to perform Bayesian inference on the parameters of a continuous-time stochastic process using CompetingClocks.jl and Turing.jl.

## Statistical Framework

### The Model

Consider a continuous-time birth-death process where, at time ``t`` with population ``N(t)``:

- **Births** occur at rate ``\beta N(t)``, i.e., waiting times are ``\text{Exponential}(1/(\beta N))``
- **Deaths** follow independent ``\text{Gamma}(k, \theta)`` lifetimes for each individual

This is a generalized semi-Markov process (GSMP) because death times are non-exponential.

### The Inference Problem

Given a fully observed event path on ``[0, T_{\text{end}}]``:

```math
\mathcal{X} = \{(e_i, t_i)\}_{i=1}^m
```

where ``e_i \in \{(\text{birth}, \cdot), (\text{death}, \text{id})\}`` and ``t_1 < t_2 < \cdots < t_m``, we want the posterior:

```math
p(\beta, k, \theta \mid \mathcal{X}) \propto p(\mathcal{X} \mid \beta, k, \theta) \, p(\beta, k, \theta)
```

### Path Likelihood

The path likelihood for a GSMP is:

```math
p(\mathcal{X} \mid \theta) = \prod_{i=1}^{m} h_{e_i}(t_i) \cdot \prod_{j \in \text{enabled}} S_j(t_i^- \to t_i) \cdot \prod_{j \in \text{enabled at } t_m} S_j(t_m \to T_{\text{end}})
```

where ``h_j(t)`` is the hazard of clock ``j`` and ``S_j(a \to b)`` is the survival probability over ``[a,b]``.

CompetingClocks computes this via `pathloglikelihood(sampler, T_end)` when the sampler is constructed with `path_likelihood=true`.

## Integration Pattern

Turing requires a `logpdf` function for likelihood terms. We wrap the CompetingClocks path likelihood in a custom `Distributions.jl` type:

```julia
struct BirthDeathPathDist <: ContinuousMultivariateDistribution
    β::Float64
    shape::Float64
    scale::Float64
    N0::Int
end
```

The `logpdf` method:
1. Creates a `SamplingContext` with `path_likelihood=true`
2. Initializes the model state and clocks
3. Replays each observed event via `fire!` and state updates
4. Returns `pathloglikelihood(sampler, T_end)`

From Turing's perspective, the observed path is a single draw from this distribution:

```julia
path ~ BirthDeathPathDist(β, shape, scale, N0)
```

## Complete Example

The full working example is available at [`examples/turing_dist.jl`](https://github.com/adolgert/CompetingClocks.jl/blob/main/examples/turing_dist.jl). Key components are shown below.

### Model State and Events

```julia
using Random, Distributions, Turing, CompetingClocks

const ClockKey = Tuple{Symbol,Int}

struct BDParams
    β::Float64       # per-individual birth rate
    shape::Float64   # Gamma shape for lifetime
    scale::Float64   # Gamma scale for lifetime
end

mutable struct BDState
    population::Set{Int}
    next_id::Int
end

BDState(N0::Int) = BDState(Set(1:N0), N0 + 1)
```

### Clock Management

```julia
function initialize!(state::BDState, sampler, params::BDParams)
    N = length(state.population)
    N == 0 && return
    enable!(sampler, (:birth, 0), Exponential(1 / (params.β * N)))
    for id in state.population
        enable!(sampler, (:death, id), Gamma(params.shape, params.scale))
    end
end

function handle_event!(state::BDState, sampler, evt::ClockKey, t::Float64, params::BDParams)
    kind, idx = evt
    if kind === :birth
        new_id = state.next_id
        state.next_id += 1
        push!(state.population, new_id)
        enable!(sampler, (:death, new_id), Gamma(params.shape, params.scale))
    elseif kind === :death
        delete!(state.population, idx)
    end

    if !isempty(state.population)
        enable!(sampler, (:birth, 0), Exponential(1 / (params.β * length(state.population))))
    else
        disable!(sampler, (:birth, 0))
    end
end
```

### Path Representation

```julia
struct BDPath
    events::Vector{ClockKey}
    times::Vector{Float64}
    T_end::Float64
end

Base.length(p::BDPath) = length(p.events)
Base.isnan(p::BDPath) = any(isnan, p.times) || isnan(p.T_end)  # Required by Turing
```

### Simulation

```julia
function simulate_path(params::BDParams; N0::Int=10, T_end::Float64=10.0, rng=Random.Xoshiro(42))
    state = BDState(N0)
    sampler = SamplingContext(ClockKey, Float64, rng)
    initialize!(state, sampler, params)

    events, times = ClockKey[], Float64[]
    while true
        when, which = next(sampler)
        (isnothing(which) || when > T_end) && break
        push!(events, which); push!(times, when)
        fire!(sampler, which, when)
        handle_event!(state, sampler, which, when, params)
    end
    BDPath(events, times, T_end)
end
```

### Custom Distribution for Turing

```julia
struct BirthDeathPathDist <: Distributions.ContinuousMultivariateDistribution
    β::Float64
    shape::Float64
    scale::Float64
    N0::Int
end

Distributions.length(d::BirthDeathPathDist) = 1

function Distributions.logpdf(d::BirthDeathPathDist, path::BDPath)
    params = BDParams(d.β, d.shape, d.scale)
    rng = Random.Xoshiro(0)  # RNG unused since we only replay
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)

    state = BDState(d.N0)
    initialize!(state, sampler, params)

    for (evt, t) in zip(path.events, path.times)
        fire!(sampler, evt, t)
        handle_event!(state, sampler, evt, t, params)
    end

    ll = pathloglikelihood(sampler, path.T_end)
    return ll isa AbstractVector ? ll[1] : ll
end

# Required by Turing for likelihood accumulation
Distributions.loglikelihood(d::BirthDeathPathDist, path::BDPath) = Distributions.logpdf(d, path)
```

### Turing Model

```julia
@model function birthdeath_model(path::BDPath, N0::Int)
    β     ~ LogNormal(log(0.3), 0.5)
    shape ~ LogNormal(log(2.0), 0.5)
    scale ~ LogNormal(log(0.5), 0.5)
    path  ~ BirthDeathPathDist(β, shape, scale, N0)
end
```

### Running Inference

```julia
# True parameters: ensure birth rate < death rate for finite simulation
# With Gamma(2, 0.5), mean lifetime = 1.0, so death rate ≈ 1.0
# Birth rate β = 0.3 < 1.0 gives declining population
true_params = BDParams(0.3, 2.0, 0.5)
observed_path = simulate_path(true_params; N0=20, T_end=10.0, rng=Random.Xoshiro(1234))

model = birthdeath_model(observed_path, 20)
chain = sample(model, MH(), 500; rng=Random.Xoshiro(2024))
```

**Example output:**
```
Parameter recovery:
  β: true=0.300, estimated=0.263 ± 0.095
  shape: true=2.000, estimated=2.073 ± 0.428
  scale: true=0.500, estimated=0.545 ± 0.162
```

## Alternative: `@addlogprob!`

Instead of defining a custom distribution, you can add the log-likelihood directly:

```julia
@model function birthdeath_model_addlogprob(path::BDPath, N0::Int)
    β     ~ LogNormal(log(0.3), 0.5)
    shape ~ LogNormal(log(2.0), 0.5)
    scale ~ LogNormal(log(0.5), 0.5)

    params = BDParams(β, shape, scale)
    sampler = SamplingContext(ClockKey, Float64, Random.Xoshiro(0); path_likelihood=true)
    state = BDState(N0)
    initialize!(state, sampler, params)

    for (evt, t) in zip(path.events, path.times)
        fire!(sampler, evt, t)
        handle_event!(state, sampler, evt, t, params)
    end

    Turing.@addlogprob! pathloglikelihood(sampler, path.T_end)
end
```

## Notes

**Parameter stability:** For birth-death processes, ensure `β < 1/(shape × scale)` to avoid population explosion during simulation. Here, with `β=0.3` and `Gamma(2, 0.5)` (mean lifetime 1.0), the death rate exceeds the birth rate.

**Sampler choice:** The Metropolis-Hastings sampler (`MH()`) works well for this likelihood. Gradient-based samplers like `NUTS` require differentiable likelihoods, which may not be straightforward for path-based models.

**Multiple distributions:** CompetingClocks supports computing likelihoods under multiple parameter settings simultaneously via `likelihood_cnt > 1` in the `SamplingContext`, useful for importance sampling schemes.
