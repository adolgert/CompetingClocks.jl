# CompetingClocks as a Gen Distribution

This document shows how to treat an entire CompetingClocks trajectory as a single random choice in [Gen.jl](https://www.gen.dev/), enabling Bayesian inference over continuous-time discrete-event systems.

## Statistical Framework

### Path Likelihood for Generalized Semi-Markov Processes

Consider a continuous-time system over horizon ``[0, T]`` where:

- At time ``t``, a set of clocks ``\mathcal{K}_t`` are enabled
- Each clock ``k`` has a survival function ``S_k(\tau) = 1 - F_k(\tau)`` and hazard ``h_k(\tau) = f_k(\tau)/S_k(\tau)``
- The system produces a path ``x = \{(t_1, e_1), \ldots, (t_n, e_n)\}`` with ``0 < t_1 < \cdots < t_n \le T``

The path density with respect to Lebesgue measure on event times follows the standard GSMP form:

```math
\log p(x \mid \theta) = \sum_{i=1}^n \log h_{e_i}(t_i \mid \mathcal{H}_{t_i^-}, \theta) - \int_0^T \sum_{k \in \mathcal{K}_s} h_k(s \mid \mathcal{H}_{s^-}, \theta) \, ds
```

The first term sums log-hazards at each firing time. The second term—the integrated hazard over all enabled clocks—accounts for the probability of *not* firing before each event.

### CompetingClocks Handles the Integrals

When you construct a `SamplingContext` with `path_likelihood=true`, CompetingClocks tracks both terms automatically. After simulation:

```julia
log_prob = pathloglikelihood(sampler, T)
```

returns the exact log path likelihood, including the probability of no further events until time ``T``.

### Integration with Gen

The strategy is:

1. **Simulate** paths using CompetingClocks' `next`/`fire!`/`enable!` loop
2. **Evaluate** ``\log p(x \mid \theta)`` via `pathloglikelihood`
3. **Wrap** these as `Gen.random` and `Gen.logpdf` for a custom distribution

From Gen's perspective, the entire trajectory becomes a single continuous random choice whose density is delegated to CompetingClocks.

---

## Example: Density-Dependent Birth-Death Process

A simple linear birth-death process either explodes or dies out. For stable dynamics suitable for inference, we use a logistic birth rate:

```math
\lambda(N) = \lambda_0 \cdot N \cdot \max\left(0, 1 - \frac{N}{K}\right)
```

where ``K`` is the carrying capacity. This ensures the population fluctuates around ``K`` rather than diverging.

The complete working example is in `examples/gen_distribution.jl`.

### Basic Types

```julia
using Random
using Distributions
using CompetingClocks

# Clock keys: (:birth, 0) for births, (:death, i) for individual i's death
const ClockKey = Tuple{Symbol,Int}

# A single event in a trajectory
struct BDEvent
    time::Float64
    key::ClockKey
end

const EventPath = Vector{BDEvent}

# Population state
mutable struct BDState
    population::Set{Int}
    next_id::Int
end

BDState(N0::Int) = BDState(Set(1:N0), N0 + 1)
```

### Simulation with Path Likelihood

```julia
function simulate_bd(
    rng::AbstractRNG,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)
    # Enable path likelihood tracking
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)
    state = BDState(N0)
    path = EventPath()

    # Density-dependent birth rate
    function birth_rate(N::Int)
        N <= 0 && return 0.0
        return λ_birth * N * max(0.0, 1.0 - N / K)
    end

    # Initialize clocks
    rate = birth_rate(length(state.population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in state.population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Simulation loop
    when, which = next(sampler)
    while !isnothing(which) && when <= t_max
        fire!(sampler, which, when)
        push!(path, BDEvent(when, which))

        # Update state
        if which[1] == :birth
            new_id = state.next_id
            state.next_id += 1
            push!(state.population, new_id)
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif which[1] == :death
            delete!(state.population, which[2])
        end

        # Refresh birth clock with updated rate
        rate = birth_rate(length(state.population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end

        when, which = next(sampler)
    end

    return path, sampler
end
```

### Replay for Log-Likelihood Evaluation

To compute ``\log p(x \mid \theta)`` for an arbitrary path ``x``, we "replay" it through a fresh sampler:

```julia
function bd_path_logpdf(
    path::EventPath,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::Float64
    # Validate path structure
    if any(ev -> ev.time < 0 || ev.time > t_max, path)
        return -Inf
    end
    for i in 2:length(path)
        if path[i].time <= path[i-1].time
            return -Inf
        end
    end

    # Fresh sampler for replay (RNG unused)
    rng = Xoshiro(0)
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)
    state = BDState(N0)

    # Same birth rate function
    function birth_rate(N::Int)
        N <= 0 && return 0.0
        return λ_birth * N * max(0.0, 1.0 - N / K)
    end

    # Initialize
    rate = birth_rate(length(state.population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in state.population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Replay events
    for ev in path
        if !isenabled(sampler, ev.key)
            return -Inf  # Invalid: firing a disabled clock
        end

        fire!(sampler, ev.key, ev.time)

        if ev.key[1] == :birth
            new_id = state.next_id
            state.next_id += 1
            push!(state.population, new_id)
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif ev.key[1] == :death
            delete!(state.population, ev.key[2])
        end

        rate = birth_rate(length(state.population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end
    end

    return pathloglikelihood(sampler, t_max)
end
```

---

## Gen Distribution Wrapper

### Distribution Type

```julia
using Gen

struct BDPathDist <: Gen.Distribution{EventPath} end
const bd_path_dist = BDPathDist()
```

### Required Methods

```julia
function Gen.random(
    ::BDPathDist,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::EventPath
    rng = Random.default_rng()
    path, _ = simulate_bd(rng, t_max, λ_birth, K, death_shape, death_scale, N0)
    return path
end

function Gen.logpdf(
    ::BDPathDist,
    path::EventPath,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::Float64
    return bd_path_logpdf(path, t_max, λ_birth, K, death_shape, death_scale, N0)
end

Gen.is_discrete(::BDPathDist) = false
Gen.has_output_grad(::BDPathDist) = false
Gen.has_argument_grads(::BDPathDist) = (false, false, false, false, false, false)
```

This is sufficient for importance sampling, SMC, and Metropolis-Hastings. For gradient-based inference (HMC/NUTS), implement `logpdf_grad`.

---

## Using the Distribution in a Gen Model

### Generative Model

```julia
@gen function bd_model(t_max::Float64, K::Float64, N0::Int)
    # Priors on rate parameters
    λ_birth = @trace(gamma(2.0, 1.0), :λ_birth)        # mean = 2
    mean_death = @trace(gamma(2.0, 2.0), :mean_death)  # mean = 4

    death_shape = 2.0
    death_scale = mean_death / death_shape

    # Entire trajectory as one random choice
    path = @trace(
        bd_path_dist(t_max, λ_birth, K, death_shape, death_scale, N0),
        :path,
    )

    return (λ_birth=λ_birth, mean_death=mean_death, path=path)
end
```

### Forward Simulation

```julia
t_max = 10.0
K = 50.0
N0 = 10

trace = Gen.simulate(bd_model, (t_max, K, N0))
retval = Gen.get_retval(trace)

println("λ_birth = ", retval.λ_birth)
println("mean_death = ", retval.mean_death)
println("Events = ", length(retval.path))
println("Log-probability = ", Gen.get_score(trace))
```

### Conditioning on Observed Data

Given an observed trajectory `obs_path`:

```julia
obs = Gen.choicemap((:path, obs_path))
trace, logw = Gen.generate(bd_model, (t_max, K, N0), obs)

# trace[:λ_birth] and trace[:mean_death] are sampled from the prior
# logw is the importance weight
```

For proper posterior inference, use Gen's inference library:

```julia
# Importance sampling
traces, weights, _ = Gen.importance_sampling(bd_model, (t_max, K, N0), obs, n_samples)

# MCMC
trace, = Gen.importance_resampling(bd_model, (t_max, K, N0), obs, 100)
for i in 1:1000
    trace, = Gen.mh(trace, Gen.select(:λ_birth, :mean_death))
end
```

---

## Generalizing to Other Models

The pattern applies to any CompetingClocks simulation:

1. **Define your simulation** using `SamplingContext(...; path_likelihood=true)`
2. **Record the path** as a sequence of (time, event) pairs
3. **Write a replay function** that:
   - Creates a fresh sampler with `path_likelihood=true`
   - Replays each event with `fire!(sampler, key, time)`
   - Returns `pathloglikelihood(sampler, T)`
4. **Wrap as `Gen.Distribution`** with `random` and `logpdf`

The simulation logic (state transitions, clock distributions) is model-specific. The Gen integration remains identical.

### Examples to Adapt

- **SIR epidemics**: Infection and recovery clocks with population-dependent rates
- **Queueing systems**: Arrival and service clocks with queue-length feedback
- **Reliability models**: Component failure and repair with dependent hazards
- **Chemical kinetics**: Reaction clocks with mass-action propensities

In each case, the path likelihood integral is handled automatically by CompetingClocks.
