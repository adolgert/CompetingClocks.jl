using Random
using Distributions
using Fleck


mutable struct ConstantBirth{T <: ContinuousUnivariateDistribution}
    birth_rate::Float64
    death_distribution::T
    next_name::Int64
    alive::Int64
end


function initialize_model!(model, sampler, rng)
    enable!(sampler, 1, Exponential(1.0 / model.birth_rate), 0.0, 0.0, rng)

    initial_population = model.birth_rate * mean(model.death_distribution)

    for name_id in 1:Int(round(initial_population))
        enable!(sampler, name_id, model.death_distribution, 0.0, 0.0, rng)
        model.next_name = name_id + 1
        model.alive += 1
    end
end


function step!(model::ConstantBirth, sampler, when::Float64, which::Int, rng)
    if which == 1
        # It's a birth. Disable and enable the birth process.
        disable!(sampler, 1, when)
        enable!(sampler, 1, Exponential(1.0 / model.birth_rate), when, when, rng)

        name_id = model.next_name
        enable!(sampler, name_id, model.death_distribution, when, when, rng)
        model.next_name += 1
        model.alive += 1
    else
        # It's a death
        disable!(sampler, which, when)
        model.alive -= 1
    end
    model.alive
end


"""
The birth-death model with constant birth and non-constant rate death is
the most basic demography problem. Let's run this and ask whether it
will stay at steady state if it started at steady state.
"""
function run_constant_birth(rng, max_step = 10000)
    birth_rate = 117.0
    death_rate = Weibull(2.0, 80)
    model = ConstantBirth(birth_rate, death_rate, 2, 0)

    sampler = FirstToFire{Int}()
    initialize_model!(model, sampler, rng)
    # Begin by dropping a few events to account for warm-up.
    when = 0.0
    (when, which) = next(sampler, when, rng)
    while when < 1e4
        step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end

    # Then collect statistics.
    total = 0
    for step_idx in 1:max_step
        total += step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, when, rng)
    end
    steady_state = model.birth_rate * mean(model.death_distribution)
    observed_state = total / max_step
    (steady_state, observed_state)
end


"""
This function below will check that, as we increase the data collected,
it gets closer to the expected average with smaller standard deviation.

julia> multiple_runs(20, 100)
(8295.084022237816, 4.797754442480227e-5, 0.014832523359510309)

julia> multiple_runs(20, 1000)
(8295.084022237816, 0.0023953678719728436, 0.011503580263544711)

julia> multiple_runs(20, 10000)
(8295.084022237816, 0.000672396777082785, 0.007723221824191141)

julia> multiple_runs(20, 100000)
(8295.084022237816, 0.0002737718817670307, 0.003979880392829035)

julia> multiple_runs(20, 1000000)
(8295.084022237816, -0.0004463841629437419, 0.0018613416183046558)

julia> multiple_runs(20, 10000000)
(8295.084022237816, 0.0003425138641836902, 0.0004266251274355263)
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
