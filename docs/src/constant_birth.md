```@meta
EditURL = "constant_birth.jl"
```

# Delayed Birth-Death Process

````julia
using Random
using Distributions
using Fleck
````

[Birth-death processes](https://en.wikipedia.org/wiki/Birth%E2%80%93death_process) are a fundamental
type of stochastic process, and are the building block of many more complicated models. Here
we demonstrate how to use Fleck to build a very simple simulation of a birth-death process
where birth occurs according to an exponential (Markov) clock, but death occurs according
to a Weibull distribution. We compare ensemble results to the known stationary distribution.
Such models have been considered many times in the literature, but a recent reference is in
["Stochastic description of delayed systems"](https://doi.org/10.1098/rsta.2012.0458) by Lafuerza and Toral.

## Model structure

The model will be stored in a struct with a type parameter that is a subtype of `ContinuousUnivariateDistribution`,
which is the distribution type for the clock associated to death.
We also define a function `initialize_model!` which enables a single clock
for the birth event, and for each individual in the initial population,
enables a death clock for that individual.

Known results tell us that the mean of the stationary distribution will be
the birth rate multiplied by the average duration alive, which we use
to choose the number of individuals in the initial population. We expect
the model to fluctuate randomly around this value.

````julia
mutable struct ConstantBirth{T <: ContinuousUnivariateDistribution}
    birth_rate::Float64
    death_distribution::T
    next_name::Int64
    alive::Int64
    when::Float64
end


function initialize_model!(model, sampler, rng)
    enable!(sampler, 1, Exponential(1.0 / model.birth_rate), 0.0, 0.0, rng)

    initial_population = model.birth_rate * mean(model.death_distribution)

    for name_id in 1:Int(round(initial_population))
        past_birth = rand(rng, model.death_distribution)
        enable!(sampler, name_id, model.death_distribution, -past_birth, 0.0, rng)
        model.next_name = name_id + 1
        model.alive += 1
    end
end;
````

## Model update

There's two classes of events that can occur in this model. Birth is always
assigned to key `1`. When it fires, we disable and enable the birth process to
reset it, and then enable a death clock for the new individual. If the firing
event was death, we simply disable the clock. We return the integrated population
over time from the `step!` method to check simulation results.

````julia
function step!(model::ConstantBirth, sampler::SSA{K,T}, when::T, which::K, rng) where {K,T}
    if which == 1
        disable!(sampler, 1, when)
        enable!(sampler, 1, Exponential(1.0 / model.birth_rate), when, when, rng)

        name_id = model.next_name
        enable!(sampler, name_id, model.death_distribution, when, when, rng)
        model.next_name += 1
        model.alive += 1
    else
        disable!(sampler, which, when)
        model.alive -= 1
    end
    duration = when - model.when
    model.when = when
    model.alive * duration
end;
````

Our run function is simple. We use the `FirstToFire` sampler, but any sampler from Fleck
capable of supporting non-Exponential distributions can be used.

````julia
function run_constant_birth(rng, max_step = 10000)
    birth_rate = 117.0
    death_rate = Weibull(2.0, 80)
    model = ConstantBirth(birth_rate, death_rate, 2, 0, 0.0)

    sampler = FirstToFire{Int,Float64}()
    initialize_model!(model, sampler, rng)
    # Begin by dropping a few events to account for burn-in.
    when = 0.0
    (when, which) = next(sampler, when, rng)
    while when < 1e4
        step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end

    # Then collect statistics.
    total::Float64 = 0.0
    start_time = when
    for _ in 1:max_step
        total += step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end
    steady_state = model.birth_rate * mean(model.death_distribution)

    observed_state = total / (when - start_time)
    (steady_state, observed_state)
end;
````

## Simulation

We check below that as we increase the data collected,
it gets closer to the expected average with smaller standard deviation.

````julia
function multiple_runs(trial_cnt = 20, max_step = 1000)
    rng = Xoshiro(837100235)
    trials = zeros(Float64, trial_cnt)
    single_expected = 0.0
    Threads.@threads for trial_idx in 1:trial_cnt
        expected, observed = run_constant_birth(rng, max_step)
        trials[trial_idx] = observed
        single_expected = expected
    end
    (single_expected,
     (mean(trials) - single_expected) / single_expected,
     std(trials) / single_expected)
end

multiple_runs(20, 100)
````

````
(8295.084022237816, 0.008025814220693235, 0.018347178373091466)
````

````julia
multiple_runs(20, 1000)
````

````
(8295.084022237816, -0.0011113782851605522, 0.012126041783397188)
````

````julia
multiple_runs(20, 10000)
````

````
(8295.084022237816, -0.003383332381424693, 0.012100328263793862)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

