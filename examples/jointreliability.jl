# # Reliability of a Group
using Distributions
using Random
using Fleck
using Logging
using StatsPlots

#
# Generalized semi-Markov processes are often used for reliability modeling.
# A simple reliability model for some vehicle might say it is either
# working or in repair. There is a distribution of the working time and
# a distribution of the repair time. Some models make repairs take A
# fixed amount of time. The Fleck library can certainly model an individual,
# but let's consider a group of vehicles.
#
# We need to use ten vehicles a day, and we start with eighteen vehicles.
# As vehicles go into repair, the ten for the day are chosen from the
# remaining usable vehicles.
#
# If we think about one individual, the states are clear.
@enum IndividualState indnotset ready working broken

# For an individual, there are also some clear transitions.
const IndividualTransitions = Dict(
    :work => (ready, working),
    :done => (working, ready),
    :break => (working, broken),
    :repair => (broken, ready)
)

mutable struct Individual
    state::IndividualState
    transition_start::Float64
    ## This is how an individual remembers its total work leading to failure.
    work_age::Float64
    work_dist::LogUniform
    fail_dist::LogNormal
    repair_dist::Weibull
    Individual(work, fail, repair) = new(
        ready, 0.0, 0.0, work, fail, repair
        )
end


# It's more complicated to consider how we start multiple individuals working
# on any given day. We have to think about the change in state for the whole
# simulation. What are all combinations of starting and final states, and
# how do we ensure they reflect our use case of a motor pool?
#
# There can be unexpected corner cases, such as when an individual finishes
# work during the 6-6:15am when jobs are being scheduled. We simplify this
# by handling the two individual transitions that change behavior for the
# whole group. These are the moment that a worker fills up the work slots
# for the day and the moment a worker frees up a work slot, either by
# being done or by breaking.

mutable struct Experiment
    time::Float64
    group::Vector{Individual}
    # Each day the group tries to start `workers_max` workers.
    workers_max::Int64
    start_time::Float64
    rng::Xoshiro
    Experiment(group::Vector, crew_size::Int, rng) = new(0.0, group, crew_size, 0.01, rng)
end


function Experiment(individual_cnt::Int, crew_size::Int, rng)
    work_rate = LogUniform(.8, 0.99) # Gamma(9.0, 0.2)
    break_rate = LogNormal(1.5, 0.4)
    repair_rate = Weibull(1.0, 2.0)
    workers = [Individual(work_rate, break_rate, repair_rate) for _ in 1:individual_cnt]
    Experiment(workers, crew_size, rng)
end


key_type(::Experiment) = Tuple{Int,Symbol}
worker_cnt(experiment::Experiment) = size(experiment.group, 1)


# If every ready worker has an active `:work` transition, then there
# must be a well-defined time for that transition, even if the
# woker becomes ready between 6:00 am and 6:15 am.
function next_work_time(now, max_hour)
    hour = now - floor(now)
    if hour < max_hour
        return 0.0, max_hour - hour
    else
        return one(hour) - hour, one(hour) + max_hour - hour
    end
end


function handle_event(when, (who, transition), experiment, sampler)
    start_state, finish_state = IndividualTransitions[transition]
    individual = experiment.group[who]
    @assert individual.state == start_state
    individual.state = finish_state
    experiment.time = when
    disable!(sampler, (who, transition), when)

    ## First look at what happens to an individual.
    if start_state == working
        work_duration = when - individual.transition_start
        @debug "Adding $work_duration to $who"
        individual.work_age += work_duration
    end

    worker_cnt = count(w.state == working for w in experiment.group)
    need_workers = worker_cnt < experiment.workers_max
    max_hour = experiment.start_time

    if transition == :done
        disable!(sampler, (who, :break), when)
        if need_workers
            rate = Uniform(next_work_time(when, max_hour)...)
            enable!(sampler, (who, :work), rate, when, when, experiment.rng)
            @debug "schedule $who for $rate"
        end

    elseif transition == :repair
        if need_workers
            rate = Uniform(next_work_time(when, max_hour)...)
            enable!(sampler, (who, :work), rate, when, when, experiment.rng)
            @debug "schedule $who for $rate"
        end
    
    elseif transition == :work
        # enable :done and :break
        enable!(sampler, (who, :done), individual.work_dist, when, when, experiment.rng)
        ## Time shift this distribution to the left because it remembers
        ## the time already worked.
        past_work = when - individual.work_age
        enable!(sampler, (who, :break), individual.fail_dist, past_work, when, experiment.rng)
        @debug "schedule $who for done or break"

    elseif transition == :break
        # If you broke, you don't get to finish your work.
        disable!(sampler, (who, :done), when)
        individual.work_age = zero(Float64)
        enable!(sampler, (who, :repair), individual.repair_dist, when, when, experiment.rng)
        @debug "schedule $who for repair"

    else
        @assert finish_state ∈ (broken, working, ready)
    end
    individual.transition_start = when

    ## Then look at changes to behavior of the whole system.
    if transition == :work && worker_cnt == experiment.workers_max
        notnow = Int[]
        for too_many in [widx for (widx, w) in enumerate(experiment.group) if w.state == ready]
            ## You don't start today.
            disable!(sampler, (too_many, :work), when)
            push!(notnow, too_many)
        end
        @debug "Unscheduling $notnow"
    elseif transition ∈ (:done, :break) && worker_cnt == experiment.workers_max - 1
        rate = Uniform(next_work_time(when, max_hour)...)
        upnext = Int[]
        for next_chance in [widx for (widx, w) in enumerate(experiment.group) if w.state == ready]
            if next_chance != who
                enable!(sampler, (next_chance, :work), rate, when, when, experiment.rng)        
                push!(upnext, next_chance)
            end
        end
        @debug "scheduling $upnext for $rate"
    end
end


mutable struct Observation
    status::Array{Int64,2}
    started_today::Array{Int64,1}
    total_age::Array{Float64,1}
    broken_duration::Array{Float64,1}
    Observation(day_cnt, individual_cnt) = new(
        zeros(Int64, 2, day_cnt),
        zeros(Int64, day_cnt),
        zeros(Float64, day_cnt),
        zeros(Float64, individual_cnt)
    )
end

days(observation::Observation) = size(observation.status, 2)


function observe(experiment::Experiment, observation::Observation, when, which)
    who, transition = which
    day_idx = Int(floor(when))
    if transition == :work
        observation.started_today[day_idx + 1] += 1
    elseif transition == :repair
        observation.broken_duration[who] += when - experiment.group[who].transition_start
    end
    day_start = Int(floor(experiment.time + next_work_time(experiment.time, experiment.start_time)[1]))
    next_start = Int(floor(when + next_work_time(when, experiment.start_time)[1]))
    if day_start != next_start
        worker_cnt = count(w.state == working for w in experiment.group)
        broken_cnt = count(w.state == broken for w in experiment.group)
        work_ages = sum(w.work_age for w in experiment.group)
        for rec_idx in day_start:next_start - 1
            observation.status[1, 1 + rec_idx] = worker_cnt
            observation.status[2, 1 + rec_idx] = broken_cnt
            observation.total_age[1 + rec_idx] = work_ages
        end
    end
end


function run(experiment::Experiment, days)
    day_cnt = Int(ceil(days))
    observation = Observation(day_cnt, worker_cnt(experiment))
    sampler = FirstToFire{key_type(experiment),Float64}()
    rng = experiment.rng
    rate = Uniform(next_work_time(0.0, experiment.start_time)...)
    for initial in 1:length(experiment.group)
        enable!(sampler, (initial, :work), rate, 0.0, 0.0, rng)
    end
    when, which = next(sampler, experiment.time, rng)
    while isfinite(when) && when < days
        observe(experiment, observation, when, which)
        @debug "$when $which"
        handle_event(when, which, experiment, sampler)
        when, which = next(sampler, experiment.time, rng)
    end
    return observation
end

# Let's make an experiment and look at the distributions.

rng = Xoshiro(9378424)
experiment = Experiment(20, 10, rng)

plot(experiment.group[1].work_dist)
plot(experiment.group[1].fail_dist)
plot(experiment.group[1].repair_dist)

