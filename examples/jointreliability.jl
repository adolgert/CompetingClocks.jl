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

# The simplest model for one vehicle:
#
# * Common choices for the failure distribution are Weibull and Lognormal.
# * Repair distribution is a fixed time, so it's dirac-distributed.
#
# How do we account for how much the vehicle is used?
# We could make some kind of age function which depends on the
# number of days it is used and the amount of time each day.
# In the language of the GSMP, every day the vehicle is used,
# we will enable the failure transition. When we disable it,
# we will ensure that we remember the total length of time that
# failure transition has run so far before failing.
#

# Possible states for an individual
@enum IndividualState indnotset ready working broken

const IndividualTransitions = Dict(
    :work => (ready, working),
    :done => (working, ready),
    :break => (working, broken),
    :repair => (broken, ready)
)

mutable struct Individual
    state::IndividualState
    work_start::Float64
    work_age::Float64
end


struct Experiment
    time::Float64
    work_max::Int64
    worked_today::Int64
    group::Vector{Individual}
end


function handle_event((who, transition), when, experiment)
    start_state, finish_state = IndividualTransitions[transition]
    @assert experiment.group[who].state == start_state
    experiment.group[who].state = finish_state
    current_working_cnt = sum(1 for w in experiment.group if w.state == working)
    crossed_threshold = start_state == working && current_working_cnt == experiment.work_max - 1

    if transition == :done || transition == :broken
        experiment.group[who].work_age += when - expriment.group[who].work_start
    end

    if finish_state == ready
        need_worker = experiment.worked_today < experiment.work_max
        rate = Lognormal()
        enable!(sampler, (who, :work), rate, when, when, rng)
    else if finish_state == working
        expriment.group[who].work_start = when
        # enable :done and :break
    else if finish_state == broken
        # enable :repair
    end

    if crossed_threshold
        # enable all which are not this one and are ready for work.
    end
end


# Transitions

Failure
    enabled when ready
    disabled when working or broken
    has memory of total time ready
    rate is lognormal

Service
    enabled when working and when total(ready) < cnt
    disabled when broken or ready
    no memory
    rate is uniform between 7-7:15 each morning

StopWork
    enabled when working
    disabled in other states
    no memory
    rate is weibull

Repair
    enabled when broken
    rate is weibull
