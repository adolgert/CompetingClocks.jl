# A Vector Addition System
import Distributions
using Random: AbstractRNG

export VectorAdditionSystem, zero_state, vas_delta, vas_initial, fire!, simstep!,
    VectorAdditionModel, VectorAdditionFSM

"""
A `VectorAdditionSystem` is an older form of simulation that is a lot
like Petri nets. The state is a vector.
The system is a list of transitions. Each transition
is an array of values. Negative numbers mean the transition needs to
take this many tokens from the state, meaning the state at those indices
must be an integer at least that large. Positive numbers mean the
transition places tokens into the state. Unlike chemical simulations,
the rate need not depend on the number of combinations of species present.

The vector `rates` is a list of functions that look at the current state
and choose a rate for a transiton.
"""
struct VectorAdditionSystem{T <: Function}
    take::Array{Int, 2}  # states x transitions
    give::Array{Int, 2}  # states x transitions
    rates::Vector{T}  # length is transitions
end


mutable struct VectorAdditionModel
    # This part is immutable.
    vas::VectorAdditionSystem
    # The mutable state is below.
    state::Vector{Int}
    when::Float64
end


"""
This observes the vector addition machine. After each event, it outputs
the event and the time.
"""
function vam_event_observer(vam::VectorAdditionModel, when::Float64, event::Int)
    return (when, event)
end


"""
    VectorAdditionFSM(vam, initializer, sampler, rng)

This puts together the model and a sampler. This is what we think of as
a simulation. We will organize this like it's a finite state machine,
so it will have an initializer, a dynamics, and an observer.
"""
mutable struct VectorAdditionFSM
    # These are both state and dynamics.
    vam::VectorAdditionModel
    sampler::Any
    rng::AbstractRNG

    # The initializer, aka iota.
    initializer::Function
    is_initialized::Bool

    # Observer, aka lambda.
    observer::Function
end


function VectorAdditionFSM(vam, initializer, sampler, rng)
    VectorAdditionFSM(vam, sampler, rng, initializer, false, vam_event_observer)
end


"""
    zero_state(vas::VectorAdditionSystem)

Return a state vector of all zeros.
"""
function zero_state(vas::VectorAdditionSystem)
    zeros(Int, size(vas.take, 1))
end


"""
    vas_delta(state::Vector{Int}, vas::VectorAdditionSystem, transition_idx::Int)

Apply the state change assocaited with the transition indexed by `transition_idx`.
"""
function vas_delta(vas::VectorAdditionSystem, transition_idx)
    state_change = vas.give[:, transition_idx] - vas.take[:, transition_idx]
    let delta = state_change
        state -> begin state .+= delta end
    end
end

"""
    vas_initial(vas::VectorAdditionSystem, initial_state)

Return a function taking a single argument `state` and set it
to be the initial state vector.
"""
function vas_initial(vas::VectorAdditionSystem, initial_state)
    let initial = initial_state
        state -> begin state .= initial end
    end
end


"""
    fire!(sampler, vas::VectorAdditionSystem, state, modify_state, rng)

Fire a transition. `sampler` is a sampler for the next state.
`modify_state` is either an initialization or a time-stepping function.
For either initialization or time-stepping, this notifies the sampler about
changes to the transitions.
"""
function fire!(sampler, vas::VectorAdditionSystem, state, modify_state, now, rng)
    former = copy(state)
    modify_state(state)
    for rate_idx in eachindex(vas.rates)
        was_enabled = all(former .- vas.take[:, rate_idx] .>= 0)
        now_enabled =  all(state .- vas.take[:, rate_idx] .>= 0)
        if was_enabled && !now_enabled
            disable!(sampler, rate_idx, now)
        elseif !was_enabled && now_enabled
            enable!(sampler, rate_idx, vas.rates[rate_idx](state), now, now, rng)
        elseif was_enabled && now_enabled
            ratefunc = vas.rates[rate_idx]
            former_rate = ratefunc(former)
            current_rate = ratefunc(state)
            if former_rate != current_rate
                enable!(sampler, rate_idx, current_rate, now, now, rng)
            end  # Else don't notify because rate is the same.
        end
    end
end


"""
Tell the finite state machine to step. A finite state machine has an input
token, but the token here is always, "Go!"
"""
function simstep!(fsm::VectorAdditionFSM)
    if !fsm.is_initialized
        fire!(fsm.sampler, fsm.vam.vas, fsm.vam.state, fsm.initializer, fsm.vam.when, fsm.rng)
        fsm.is_initialized = true
    end
    (when, what) = next(fsm.sampler, fsm.vam.when, fsm.rng)
    if when < Inf
        action = vas_delta(fsm.vam.vas, what)
        fsm.vam.when = when
        fire!(fsm.sampler, fsm.vam.vas, fsm.vam.state, action, fsm.vam.when, fsm.rng)
        fsm.observer(fsm.vam, when, what)
    else
        (when, what)
    end
end
