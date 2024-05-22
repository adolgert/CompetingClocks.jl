# # Reliability Model of a Work Crew

# This is an extended example using a reliability model.

using ColorSchemes
using Distributions
using Fleck
using Logging
using Plots
using Random
using StatsPlots

# ## Overview
#
# The classic model for reliability is of a machine that is either
# working or broken. There is a distribution of failure times and a distribution
# of repair times [1]. Let's extend this idea to the reliability of a vehicle
# motor pool.
#
#  1. There are 16 vehicles.
#  2. Every morning, 10 vehicles go out for work. They all leave in the first 15 mins.
#  3. Each vehicle works at least half a day, at most the whole day.
#  4. While working, a vehicle can break, where the chance of breaking depends
#     on the total time since it was last repaired.
#  5. There is a distribution of repair times.
#
# Number two, above, says that vehicles start in the first 15 minutes. This
# plan will initiate ten transitions in rapid succession. We could,
# instead, start all ten vehicles at the same time, using a single transition.
# Either would work.
#
# Fleck will take care of the timing of all of the events, but we will see that
# there is plenty of work to track the state of all of the vehicles. This
# extended example shows that, if we wanted to create more reliability models,
# it would make sense to create a framework for reliability modeling, one that
# uses Fleck underneath.
#
# ## Define State for the Model
#
# If we think about an individual vehicle, the states are ready, working, or broken.
#
@enum IndividualState ready working broken

# There are four allowed tansitions among the three states because a vehicle
# breaks only while it's working.
#
const IndividualTransitions = Dict(
    :work => (ready, working),
    :done => (working, ready),
    :break => (working, broken),
    :repair => (broken, ready)
);
#
# An individual has state and parameters. In the language of generalized
# semi-Markov processes, this state is called the *physical state* in order
# to distinguish it from the state of each enabled transition for each
# vehicle.
#
mutable struct Individual
    ## State for the individual
    state::IndividualState
    work_age::Float64 ## How an individual remembers its total work leading to breaks.
    transition_start::Float64  ## This is bookkeeping.
    ## Parameters for the individual
    done_dist::LogUniform
    fail_dist::LogNormal
    repair_dist::Weibull
    Individual(work, fail, repair) = new(
        ready, 0.0, 0.0, work, fail, repair
        )
end
#
# The simulation as a whole is the state of the individuals and the system time.
# We put some parameters here:
#
#  * `workers_max` - Each morning, this many vehicles start driving, if at least
#    this many vehicles are ready, instead of broken.
#  * `start_time` - Vehicles start in the first 15 min or so, and this is that 15 min.
#
mutable struct Experiment
    time::Float64
    group::Vector{Individual}
    ## Each day the group tries to start `workers_max` workers.
    workers_max::Int64
    start_time::Float64
    rng::Xoshiro
    Experiment(group::Vector, crew_size::Int, rng) = new(0.0, group, crew_size, 0.01, rng)
end
#
# Make a simulation by making individuals.
#
function Experiment(individual_cnt::Int, crew_size::Int, rng)
    done_rate = LogUniform(.8, 0.99) # Gamma(9.0, 0.2)
    break_rate = LogNormal(1.5, 0.4)
    repair_rate = Weibull(1.0, 2.0)
    workers = [Individual(done_rate, break_rate, repair_rate) for _ in 1:individual_cnt]
    Experiment(workers, crew_size, rng)
end
#
# And make some helpers. The `key_type` says that we will track transitions using
# a tuple of (index of vehicle, symbol to identify the transition).
#
key_type(::Experiment) = Tuple{Int,Symbol};
worker_cnt(experiment::Experiment) = size(experiment.group, 1);
#
# ## Define Transitions for the Model
#
# If we were modeling one individual, transitions would be very simple,
# but by asking that ten vehicles work every morning, we require that those
# individuals interact.
#
# One way to think clearly about interactions is to think about the state
# of the whole system. If less than ten vehicles are currently working, then
# every ready vehicle must have an enabled transition to start work at the
# next available time. Once the tenth vehicle begins working, all of those
# transitions need to be disabled.
#
# It's implied that the start of each day happens at 1.0, 2.0, 3.0, etc. When a
# vehicle becomes ready, or when the total working vehicles drops below
# ten, then each ready vehicle could work at a future time. This function
# takes in the current time and returns two times, relative to the current
# time, between which the vehicle can start work.
#
function next_work_time(now, fifteen_minutes)
    hour = now - floor(now)
    if hour < fifteen_minutes ## If vehicles are still going out today.
        return 0.0, fifteen_minutes - hour
    else ## You can't start until tomorrow.
        return one(hour) - hour, one(hour) + fifteen_minutes - hour
    end
end;
#
# Now we handle simulation events. This function's complexity is an argument
# for using a framework like a queueing model, a generalized stochastic Petri
# net, or some other continuous-time simulation framework.
#
# The arguments are:
# 
#  * `when` - The time of the next event.
#  * `(who, transition)` - This expands the `key_type`, which identifies the transition.
#  * `experiment` - It's our simulation data.
#  * `sampler` - This is a [Fleck.SSA](@ref) from Fleck to enable and disable transitions.
#
# The first few statements of the function are automatic for any transition.
# Then this handler works through the transition types.
function handle_event(when, (who, transition), experiment, sampler)
    start_state, finish_state = IndividualTransitions[transition]
    individual = experiment.group[who]
    @assert individual.state == start_state
    individual.state = finish_state
    experiment.time = when
    disable!(sampler, (who, transition), when)

    ## If a vehicle is done work, or if they break, then include the time worked
    ## in their total work age.
    if start_state == working
        work_duration = when - individual.transition_start
        @debug "Adding $work_duration to $who"
        individual.work_age += work_duration
    end

    ## The state of the system, as a whole, depends on the total number
    ## currently working.
    worker_cnt = count(w.state == working for w in experiment.group)
    need_workers = worker_cnt < experiment.workers_max
    max_hour = experiment.start_time

    ## When an individual was working, there were two possible transitions,
    ## one to `ready`, and one to `broken`. Don't forget to disable the `:break`
    ## transition. Then schedule the next day's work only if the system has less
    ## than ten working.
    if transition == :done
        disable!(sampler, (who, :break), when)
        if need_workers
            rate = Uniform(next_work_time(when, max_hour)...)
            enable!(sampler, (who, :work), rate, when, when, experiment.rng)
            @debug "schedule $who for $rate"
        end

    ## A `:repair` transition can happen at any time, including during the first
    ## fifteen minutes of a day.
    elseif transition == :repair
        if need_workers
            rate = Uniform(next_work_time(when, max_hour)...)
            enable!(sampler, (who, :work), rate, when, when, experiment.rng)
            @debug "schedule $who for $rate"
        end

    ## The `:work` transition represents a vehicle going out to work for the day.
    ## This enables two possible transitions, finishing work or breaking. The
    ## breaking transition is interesting because it has what Zimmerman [2] calls
    ## "memory." It remembers how long it was previously enabled
    elseif transition == :work
        ## enable :done and :break
        enable!(sampler, (who, :done), individual.done_dist, when, when, experiment.rng)
        ## Time shift this distribution to the left because it remembers
        ## the time already worked.
        past_work = when - individual.work_age
        enable!(sampler, (who, :break), individual.fail_dist, past_work, when, experiment.rng)
        @debug "schedule $who for done or break"

    ## When a vehicle breaks, the only option is to repair it. This resets the work age.
    elseif transition == :break
        ## If you broke, you don't get to finish your work.
        disable!(sampler, (who, :done), when)
        individual.work_age = zero(Float64)
        enable!(sampler, (who, :repair), individual.repair_dist, when, when, experiment.rng)
        @debug "schedule $who for repair"

    else
        @assert transition ∈ keys(IndividualTransitions)
    end
    individual.transition_start = when

    ## We haven't handled how we ensure that at most ten vehicles start work every
    ## morning. For that, we need to think about the system as a whole, explicitly
    ## by looking at the current worker count and whether it crossed the threshold of
    ## ten workers.
    ##
    ## If a vehicle just started and is the tenth worker, then cancel the ability of
    ## all other vehicles to work.
    if transition == :work && worker_cnt == experiment.workers_max
        notnow = Int[]
        for too_many in [widx for (widx, w) in enumerate(experiment.group) if w.state == ready]
            ## You don't start today.
            disable!(sampler, (too_many, :work), when)
            push!(notnow, too_many)
        end
        @debug "Unscheduling $notnow"

    ## If a vehicle stopped work, either by finishing or breaking, and it was the
    ## first of the work crew to quit, then notify all `ready` vehicles that they
    ## should start work at the start of the next morning.
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
end;
#
# ## Configure the Model
#
# For anything other than an example, the most important step would be
# configuring the model so that it matches observations. Here, however, we
# have put this directly into the `Experiment` type. Here is a plot of
# the distributions.
#
function show_distributions()
    experiment = Experiment(16, 10, Xoshiro(9378424))

    plot(experiment.group[1].done_dist, label="Done")
    plot!(experiment.group[1].fail_dist, label="Break")
    plot!(experiment.group[1].repair_dist, label="Repair")
    title!("Distributions for Transitions")
end
show_distributions()
# The short blue line in the upper-left is the probability distribution
# function for a `LogUniform` distribution that represents the time
# a vehicle drives on a single day. You can see that these vehicles
# break about once a week and take a couple of days to repair, on average.

# ## Run the Simulation
#
# Running a simulation means we are sampling from the stochastic process.
# For a continuous-time stochastic process like this, that means asking the
# sampler when the next transition is and which transition it is. If there
# are no possible transitions, the next time will be infinite and the chosen
# transition will be `nothing`.
#
function run(experiment::Experiment, observation, days)
    sampler = FirstToFire{key_type(experiment),Float64}()
    rng = experiment.rng
    rate = Uniform(next_work_time(0.0, experiment.start_time)...)
    for initial in 1:length(experiment.group)
        enable!(sampler, (initial, :work), rate, 0.0, 0.0, rng)
    end
    when, which = next(sampler, experiment.time, rng)
    while isfinite(when) && when < days
        ## We use different observers to record the simulation.
        observe(experiment, observation, when, which)
        @debug "$when $which"
        handle_event(when, which, experiment, sampler)
        when, which = next(sampler, experiment.time, rng)
    end
end
#
# ## Observers
#
# Without care, data collection from continuous-time simulation can generate a lot of data quickly.
# In many cases, especially for performance analysis, not every event time and transition 
# is important. Therefore, to avoid saving the raw data stream, we use observers of the system to summarize that data.
# Construction of observers is important as it connects the simulation to tools or analyses
# that may be used to guide decision-making. It also identifies which variables and metrics
# are most important.
# 
# Let's look at a few examples for the vehicle crew.
#
# ### Continuous-time Summary Observer
#
# The first example will record data at every transition, but it will record
# only the total number of working or broken vehicles.
#
# This represents a single time point.
struct ContinuousRec
    working::Int64
    broken::Int64
    total_age::Float64
    time::Float64
end

# The observer stores a vector of those single time points.
mutable struct ObserveContinuous
    state::Vector{ContinuousRec}
    ObserveContinuous() = new([ContinuousRec(0, 0, 0.0, 0.0)])
end
#
# This observer keeps a running sum of the number working and broken.
# Note that it has to know how different transitions change those numbers.
# The relationship between the transition and how it changes counts is called
# stochiometry (or stoichiometry), because it was first observed for chemical simulations.
# Both chemical simualtions and GSPN would have this information encoded in
# a formal model.
#
function observe(experiment::Experiment, observation::ObserveContinuous, when, which)
    who, transition = which
    working = observation.state[end].working
    broken = observation.state[end].broken
    if transition == :work
        working += 1
    elseif transition == :done
        working -= 1
    elseif transition == :break
        broken += 1
        working -= 1
    elseif transition == :repair
        broken -= 1
    else
        @assert transition ∈ (:work, :done, :break, :repair)
    end
    total_age = sum(w.work_age for w in experiment.group)
    push!(observation.state, ContinuousRec(working, broken, total_age, when))
end;
#
# This plot shows a timeline of the count of working and broken vehicles
# over five days.
#
function plot_timeline(obs::ObserveContinuous, experiment::Experiment)
    state = obs.state
    last_time = state[end].time
    first_idx = findlast([x.time < last_time - 5 for x in state])
    times = [x.time for x in state[first_idx:end]]
    times .-= 4009
    working = [x.working for x in state[first_idx:end]]
    broken = [x.broken for x in state[first_idx:end]]
    ready = [worker_cnt(experiment) - x.working - x.broken for x in state[first_idx:end]]
    plot(times, working, label="working", line=(:steppost, 2))
    plot!(times, broken, label="broken", line=(:steppost, 2))
    xlabel!("Time [days]")
    ylabel!("Status [count]")
    title!("Timeline of Work Crew Over Five Days")
end

function show_typical_timeline()
    rng = Xoshiro(9234232)
    years = 11
    day_cnt = 365 * years
    worker_cnt = 16
    experiment = Experiment(worker_cnt, 10, rng)
    observation = ObserveContinuous()
    run(experiment, observation, day_cnt)

    plot_timeline(observation, experiment)
end
show_typical_timeline()
#
#
# ### Once-a-day Observation of the Working State
#
# Suppose that we want to match our simulation to observation data that counts,
# every day, how many vehicles went out and how many were broken that morning.
# This observer records that status each day.
#
mutable struct ObserveLots
    status::Array{Int64,2}
    started_today::Array{Int64,1}
    total_age::Array{Float64,1}
    broken_duration::Array{Float64,1}
    ObserveLots(day_cnt, individual_cnt) = new(
        zeros(Int64, 2, day_cnt),
        zeros(Int64, day_cnt),
        zeros(Float64, day_cnt),
        zeros(Float64, individual_cnt)
    )
end

days(observation::ObserveLots) = size(observation.status, 2);
#
# This observer waits until the current transition time is just after the
# first 15min of the day. Then it records every vehicle's status.
#
function observe(experiment::Experiment, observation::ObserveLots, when, which)
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
end;
#
# Now we can use this observer to make a plot.
#
function walk_simulation()
    day_cnt = 20
    experiment = Experiment(16, 10, Xoshiro(979798))
    observation = ObserveLots(day_cnt, worker_cnt(experiment))
    run(experiment, observation, day_cnt)
    plot(1:day_cnt, observation.status[2, :], seriestype=:scatter, label="repair",
        yticks=0:2:10)
    plot!(1:day_cnt, observation.status[1, :], seriestype=:scatter, label="working")
    title!("Number Working or in Repair")
end
walk_simulation()
#
# ### Distribution of Broken Vehicles and Probability of Missing Crew
#
# If we were focused more on the small probability that there wouldn't be
# enough vehicles in the morning to start a full ten, then we want to
# understand the histogram of how many vehicles are broken on any given
# day.
#
mutable struct ObserveHistogram
    counts::Array{Int64,2}
    working::Int64
    broken::Int64
    burn::Float64
    ObserveHistogram(e::Experiment, burn) = new(
        zeros(Int64, e.workers_max + 1, worker_cnt(e) + 1), 0, 0, burn)
end
must_work(o::ObserveHistogram) = size(o.counts, 1)
total_workers(o::ObserveHistogram) = size(o.counts, 2)

function observe(experiment::Experiment, observation::ObserveHistogram, when, which)
    if when > observation.burn
        day_start = Int(floor(experiment.time + next_work_time(experiment.time, experiment.start_time)[1]))
        next_start = Int(floor(when + next_work_time(when, experiment.start_time)[1]))
        if day_start != next_start
            observation.counts[observation.working + 1, observation.broken + 1] += 1
        end
    end

    who, transition = which
    if transition == :work
        observation.working += 1
    elseif transition == :done
        observation.working -= 1
    elseif transition == :break
        observation.broken += 1
        observation.working -= 1
    elseif transition == :repair
        observation.broken -= 1
    else
        @assert transition ∈ (:work, :done, :break, :repair)
    end
end;
#
# This observer will help us see what happens if we keep the same total
# number of vehicles but send more out each day for work.
#
function compare_across_workers(obs::Vector{ObserveHistogram}, labels, title)
    firstplot = true
    cols = palette(:tableau_20, length(obs))
    for (obs_idx, observation) in enumerate(obs)
        worker_cnt = total_workers(observation)
        broken = vec(sum(observation.counts, dims=1))
        cnt = findlast(broken .> 0)
        normed = broken[1:cnt] / sum(broken)
        if firstplot
            firstplot = false
            plot(1:cnt, normed, color=cols[obs_idx], seriestype=:scatter, markersize=2.5, label=false)
            plot!(1:cnt, normed, color=cols[obs_idx], label=labels[obs_idx], legendtitle="Crew Size")
        else
            plot!(1:cnt, normed, color=cols[obs_idx], seriestype=:scatter, markersize=2.5, label=false)
            plot!(1:cnt, normed, color=cols[obs_idx], label=labels[obs_idx])
        end
    end
    xlabel!("Count of Broken")
    ylabel!("Probability Mass")
    title!(title)
end

function show_competition_effect()
    rng = Xoshiro(4377124)
    observations = ObserveHistogram[]
    labels = String[]
    years = 10
    day_cnt = 365 * years
    worker_cnt = 20
    for must_work in [1, 5, 10, 15, 20]
        experiment = Experiment(worker_cnt, must_work, rng)
        burn = min(day_cnt ÷ 10, 3650)
        observation = ObserveHistogram(experiment, burn)
        run(experiment, observation, day_cnt)
        push!(observations, observation)
        push!(labels, string(must_work))
    end

    compare_across_workers(observations, labels, "Number Broken as Crew Increases")
end
show_competition_effect()
#
# ## References
#
# 1. Limnios, Nikolaos, and Gheorghe Oprisan. Semi-Markov processes and
#    reliability. Springer Science & Business Media, 2012.
#
# 2. Zimmermann, Armin. Stochastic discrete event systems. Springer, Berlin
#    Heidelberg New York, 2007.
