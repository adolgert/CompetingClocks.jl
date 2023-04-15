using Random
using Distributions
using Fleck


"""
This will be a very simple model. One process is birth, with a constant
rate, so it's an exponential distribution. The other is death, with
a rate that is determined by a Weibull distribution. We expect the steady
state to be the rate of birth times the average duration alive, according
to the Weibull distribution.
"""
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

    # Start out with some born.
    for name_id in 1:Int(round(initial_population))
        # We initialize the current state as though it had been running and particles
        # were born in the past. This system should _start_ in steady-state.
        past_birth = rand(rng, model.death_distribution)
        enable!(sampler, name_id, model.death_distribution, -past_birth, 0.0, rng)
        model.next_name = name_id + 1
        model.alive += 1
    end
end


function step!(model::ConstantBirth, sampler, when::Float64, which::Int, rng)
    if which == 1
        # It's a birth. Disable and enable the birth process to reset it.
        disable!(sampler, 1, when)
        enable!(sampler, 1, Exponential(1.0 / model.birth_rate), when, when, rng)

        # And register the newly-born for death.
        name_id = model.next_name
        enable!(sampler, name_id, model.death_distribution, when, when, rng)
        model.next_name += 1
        model.alive += 1
    else
        # It's a death
        disable!(sampler, which, when)
        model.alive -= 1
    end
    duration = when - model.when
    model.when = when
    model.alive * duration  # Integrated population over time.
end


"""
The birth-death model with constant birth and non-constant rate death is
the most basic demography problem. Let's run this and ask whether it
will stay at steady state if it started at steady state.
"""
function run_constant_birth(rng, max_step = 10000)
    birth_rate = 117.0
    death_rate = Weibull(2.0, 80)
    model = ConstantBirth(birth_rate, death_rate, 2, 0, 0.0)

    # sampler = FirstToFire{Int}()
    sampler = NextReaction{Int}()
    initialize_model!(model, sampler, rng)
    # Begin by dropping a few events to account for warm-up.
    when = 0.0
    (when, which) = next(sampler, when, rng)
    while when < 1e4
        step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end

    # Then collect statistics.
    total::Float64 = 0.0
    start_time = when
    for step_idx in 1:max_step
        total += step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end
    steady_state = model.birth_rate * mean(model.death_distribution)
    
    observed_state = total / (when - start_time)
    (steady_state, observed_state)
end


"""
This function below will check that, as we increase the data collected,
it gets closer to the expected average with smaller standard deviation.

julia> multiple_runs(20, 100)
(8295.084022237816, 0.011846119314536163, 0.02337492577246681)

julia> multiple_runs(20, 1000)
(8295.084022237816, -0.00411412768534268, 0.012532715853043242)

julia> multiple_runs(20, 10000)
(8295.084022237816, -0.001749444764794219, 0.009215137332073057)

julia> multiple_runs(20, 100000)
(8295.084022237816, 0.0006730212251042666, 0.005737445235070591)

julia> multiple_runs(20, 1000000)
(8295.084022237816, -0.0005756412929582505, 0.0016894450641947503)

julia> multiple_runs(20, 10000000)
(8295.084022237816, 4.9351444569388426e-5, 0.0005657842092560583)

That looks nice, right?
"""
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
