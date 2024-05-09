# # Reliability of a Group
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
    work_start::Float64
    ## This is how an individual remembers its total work leading to failure.
    work_age::Float64
end

# It's more complicated to consider how we start multiple individuals working
# on any given day. We have to think about the change in state for the whole
# simulation. What are all combinations of starting and final states, and
# how do we ensure they reflect our use case of a motor pool?
#
# The basic idea is to create a transition for all ready individuals who
# could start work, and to time that transition to be in a narrow window, say
# 6-6:15 am, or between 0.25 and 0.26, if we use 1 day as the time unit.
# There can be weird corner cases. For instance, what if one individual
# starts work at 0.25 and completes work at 0.255? Does that count as work
# for the day? What happens if work lasts more than one day?
#
# We will simplify this by defining invariants. If we know what is always
# true, then it's easier to set up transitions to guarantee that.
#
# 1. Every ready worker has an active :work transition. If it isn't
#    set to start today, then it is set to start tomorrow.
# 

struct Experiment
    time::Float64
    group::Vector{Individual}
    # Each day the group tries to start `workers_max` workers.
    workers_max::Int64
    start_times::Tuple{Float64,Float64}
    rng::Xoshiro
end


function next_work_time(now, min_hour, max_hour, need_workers)
    day = floor(now)
    hour = now - day
    if hour < min_hour
        return (day + min_hour, day + max_hour)
    else if hour < max_hour
        if need_workers
            return (day + min_hour, day + max_hour)
        else
            return (day + one(day) + min_hour, day + one(day) + max_hour)
        end
    else
        return (day + one(day) + min_hour, day + one(day) + max_hour)
    end
end


function handle_event((who, transition), when, experiment, sampler)
    start_state, finish_state = IndividualTransitions[transition]
    @assert experiment.group[who].state == start_state
    experiment.group[who].state = finish_state

    if transition == :done || transition == :broken
        experiment.group[who].work_age += when - experiment.group[who].work_start
    end

    worker_cnt = count(w.state == working for w in experiment.group)
    need_workers = worker_cnt < experiment.workers_max
    min_hour, max_hour = experiment.group.start_times

    if finish_state == ready
        rate = Uniform(next_work_time(now, min_hour, max_hour, need_workers)...)
        enable!(sampler, (who, :work), rate, when, when, experiment.rng)    

    else if finish_state == working
        expriment.group[who].work_start = when
        # enable :done and :break
        work_rate = lognormal()
        enable!(sampler, (who, :done), work_rate, when, when, experiment.rng)
        break_rate = lognormal()
        ## Time shift this distribution to the left because it remembers
        ## the time already worked.
        past_work = when - experiment.group[who].work_age
        enable!(sampler, (who, :break), break_rate, when, past_work, experiment.rng)
        if !need_workers
            for off_day in [w for w in experiment.group if w.state == ready]
                ## You don't start today.
                disable!(sampler, (off_day, :work), when)
                ## You might start tomorrow.
                rate = Uniform(next_work_time(now, min_hour, max_hour, need_workers)...)
                enable!(sampler, (who, :work), rate, when, when, experiment.rng)        
            end
        end

    else if finish_state == broken
        experiment.group[who].work_age = zero(Float64)
        repair_rate = Weibull()
        enable!(sampler, (who, :repair), repair_rate, when, when, experiment.rng)
    end
end
