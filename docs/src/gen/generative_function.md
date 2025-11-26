# CompetingClocks as a Gen Generative Function

This document shows how to wrap a CompetingClocks simulator as a custom [Gen.jl](https://www.gen.dev/) generative function, enabling it to serve as a building block in larger probabilistic models where the event path is latent.

## Statistical Framework

### Path Likelihood for Black-Box Simulators

A CompetingClocks simulator defines a probability distribution over event paths:

- **Parameters**: ``\theta`` (rates, shape parameters, etc.)
- **Event path**: ``\omega = \{(t_1, e_1), \ldots, (t_n, e_n)\}`` over horizon ``[0, T]``
- **Path density**: ``p_\theta(\omega)``

With `path_likelihood=true`, CompetingClocks computes the exact log-density via `pathloglikelihood(sampler, T)`. This includes both the probability of observed events and the survival probability (no further events before ``T``).

### The Generative Function Interface

Gen's generative function interface (GFI) wraps simulators that:

1. Take arguments (parameters, initial conditions, horizon)
2. Use internal randomness to produce a trace
3. Return a value (here: a summary like final population)
4. Report a log-probability score for the trace

By implementing `simulate`, `generate`, and trace accessors, any CompetingClocks model becomes a first-class Gen building block.

### When to Use This Approach

The generative function approach is appropriate when:

- The event path is **latent**—you observe consequences (final state, counts) rather than the path itself
- You want to **embed** the simulator in a larger model with priors and observation likelihoods
- You need Gen's inference machinery (importance sampling, MCMC) over model parameters

For conditioning directly on observed paths, see [CompetingClocks as a Gen Distribution](distribution.md).

---

## Example: Density-Dependent Birth-Death Process

A linear birth-death process either explodes or dies out. For stable dynamics suitable for inference, we use a logistic birth rate:

```math
\lambda(N) = \lambda_0 \cdot N \cdot \max\left(0, 1 - \frac{N}{K}\right)
```

where ``K`` is the carrying capacity. The population fluctuates around ``K`` rather than diverging.

The complete working example is in `examples/gen_generative.jl`.

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
```

### Simulation with Path Likelihood

```julia
function simulate_bd_path(
    birth_rate::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    init_pop::Int,
    tmax::Float64,
    rng::AbstractRNG,
)
    # Initial state: individuals are numbered 1:init_pop
    population = Set(1:init_pop)
    next_id = init_pop + 1

    # CompetingClocks sampler with path likelihood enabled
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)

    # Density-dependent birth rate: logistic model
    function birth_rate_total(N::Int)
        N <= 0 && return 0.0
        return birth_rate * N * max(0.0, 1.0 - N / K)
    end

    # Enable initial birth and death clocks
    rate = birth_rate_total(length(population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Simulate until horizon tmax or extinction
    path = BDEvent[]
    when, which = next(sampler)

    while !isnothing(which) && when <= tmax && !isempty(population)
        fire!(sampler, which, when)
        push!(path, BDEvent(when, which))

        if which[1] == :birth
            new_id = next_id
            next_id += 1
            push!(population, new_id)
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif which[1] == :death
            delete!(population, which[2])
        end

        # Refresh birth clock with density-dependent rate
        rate = birth_rate_total(length(population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end

        when, which = next(sampler)
    end

    # Log-likelihood of the whole path, including survival up to tmax
    logp = pathloglikelihood(sampler, tmax)

    return path, length(population), logp
end
```

---

## Gen Generative Function Wrapper

### Trace Type

The trace stores the simulation output and implements Gen's trace interface:

```julia
using Gen

struct BDTrace <: Gen.Trace
    args::Tuple                # (birth_rate, K, death_shape, death_scale, init_pop, tmax)
    choices::Gen.ChoiceMap     # :path => Vector{BDEvent}
    retval::Int                # final population size
    logp::Float64              # log p(path | args)
    gen_fn::Any                # back-pointer to the generative function
end
```

### Generative Function Type

The generative function wraps the simulator and holds an RNG:

```julia
struct BDPathGF{R<:AbstractRNG} <: Gen.GenerativeFunction{Int,BDTrace}
    rng::R
end

BDPathGF(rng::AbstractRNG) = BDPathGF{typeof(rng)}(rng)
BDPathGF() = BDPathGF(Random.default_rng())
```

### Required GFI Methods

```julia
function Gen.simulate(gen_fn::BDPathGF, args::Tuple)
    birth_rate, K, death_shape, death_scale, init_pop, tmax = args

    path, final_pop, logp = simulate_bd_path(
        birth_rate, K, death_shape, death_scale, init_pop, tmax, gen_fn.rng
    )

    choices = Gen.choicemap((:path, path))
    return BDTrace(args, choices, final_pop, logp, gen_fn)
end

# Default proposal q = p, so generate calls simulate with weight 0
function Gen.generate(gen_fn::BDPathGF, args::Tuple, constraints::Gen.ChoiceMap)
    trace = Gen.simulate(gen_fn, args)
    return trace, 0.0
end

Gen.generate(gen_fn::BDPathGF, args::Tuple) =
    Gen.generate(gen_fn, args, Gen.choicemap())

# Trace accessors
Gen.get_args(trace::BDTrace) = trace.args
Gen.get_retval(trace::BDTrace) = trace.retval
Gen.get_choices(trace::BDTrace) = trace.choices
Gen.get_score(trace::BDTrace) = trace.logp
Gen.get_gen_fn(trace::BDTrace) = trace.gen_fn

Base.getindex(trace::BDTrace, addr) = getindex(trace.choices, addr)

# No gradients in this example
Gen.has_argument_grads(::BDPathGF) = (false, false, false, false, false, false)
Gen.accepts_output_grad(::BDPathGF) = false

# Project returns the full log-probability
Gen.project(trace::BDTrace, ::Gen.Selection) = trace.logp
```

### Standalone Usage

The generative function can be used directly:

```julia
bd_gf = BDPathGF(Xoshiro(42))

trace = Gen.simulate(bd_gf, (2.0, 50.0, 2.0, 2.0, 10, 10.0))
final_pop = Gen.get_retval(trace)  # Int: final population
path = trace[:path]                 # Vector{BDEvent}: full trajectory
logp = Gen.get_score(trace)         # Float64: log p(path | params)
```

---

## Using the Generative Function in a Gen Model

### Inference Model

Embed the generative function in a model with priors over parameters and an observation likelihood:

```julia
const global_bd_gf = BDPathGF(Xoshiro(123))

@gen function bd_inference_model(K::Float64, tmax::Float64)
    # Prior on the birth rate (lognormal via exp of normal)
    log_birth_rate = @trace(normal(0.0, 0.5), :log_birth_rate)
    birth_rate = exp(log_birth_rate)

    # Fixed death parameters and initial population
    death_shape = 2.0
    death_scale = 2.0
    init_pop = 10

    # Draw a whole birth-death path via CompetingClocks
    final_pop = @trace(
        global_bd_gf(birth_rate, K, death_shape, death_scale, init_pop, tmax),
        :simulation
    )

    # Noisy observation of the final population
    y = @trace(normal(float(final_pop), 2.0), :y)

    return (birth_rate=birth_rate, final_pop=final_pop, y=y)
end
```

### Forward Simulation

```julia
K = 50.0
tmax = 10.0

trace = Gen.simulate(bd_inference_model, (K, tmax))
retval = Gen.get_retval(trace)

println("Sampled birth_rate: ", retval.birth_rate)
println("Final population: ", retval.final_pop)
println("Model log-probability: ", Gen.get_score(trace))
```

### Importance Sampling Inference

Condition on an observed population count and infer the birth rate:

```julia
# Observed data
observed_y = 40.0

# Condition on the observation
observations = Gen.choicemap((:y, observed_y))

# Run importance sampling
n_samples = 100
traces, log_weights, lml_est = Gen.importance_sampling(
    bd_inference_model, (K, tmax), observations, n_samples
)

# Compute importance-weighted posterior mean
posterior_birth_rates = [exp(tr[:log_birth_rate]) for tr in traces]
weights = exp.(log_weights .- maximum(log_weights))
weights ./= sum(weights)

mean_birth_rate = sum(weights .* posterior_birth_rates)
```

### Accessing the Event Path

The event path is stored in the generative function's choice map. Access it via `get_submap`:

```julia
example_trace = traces[1]
choices = Gen.get_choices(example_trace)
sim_choices = Gen.get_submap(choices, :simulation)
example_path = sim_choices[:path]
```

---

## Generalizing to Other Models

The pattern applies to any CompetingClocks simulation:

1. **Define your simulation** using `SamplingContext(...; path_likelihood=true)`
2. **Return** the path, summary statistics, and `pathloglikelihood(sampler, T)`
3. **Create a trace type** storing arguments, choices, return value, and score
4. **Implement GFI methods**: `simulate`, `generate`, and accessors

### Model-Specific Components

| Component | What to Customize |
|-----------|-------------------|
| `simulate_*_path` | State transitions, clock distributions, event handling |
| Return value | Summary statistics relevant to your observations |
| Trace type | Fields for your specific outputs |

### Extensions

- **MCMC over paths**: Implement `update`/`regenerate` to propose path modifications
- **Multiple trajectories**: Call the generative function in a loop within your model
- **Gradient-based inference**: Implement `has_argument_grads` and gradient methods
- **Biased proposals**: Use `likelihood_cnt > 1` for importance sampling with mixture proposals

### Examples to Adapt

- **SIR epidemics**: Return final recovered count; observe case reports
- **Queueing systems**: Return throughput metrics; observe service times
- **Reliability models**: Return component states; observe maintenance logs
- **Chemical kinetics**: Return species concentrations; observe assay measurements
