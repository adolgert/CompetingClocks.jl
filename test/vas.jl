# A Vector Addition System
module VectorAddition

import Distributions
using Random: AbstractRNG
using CompetingClocks

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
struct VectorAdditionSystem{Ar <: AbstractArray{<:Int,2}, Ra <: AbstractArray{<:Function}}
    take::Ar  # states x transitions
    give::Ar  # states x transitions
    rates::Ra  # length is transitions
end


# Just the state for the vector addition.
# This is mutable, while the vector addition system is immutable.
mutable struct VectorAdditionState{St <: AbstractArray{<:Int}}
    state::St  # This is X, sometimes called "physical state."
    when::Float64       # This is T
end


"""
This observes the vector addition machine. After each event, it outputs
the event and the time.
"""
function vam_event_observer(Q::VectorAdditionState, when::Float64, event::Int)
    return (when, event)
end


"""
    VectorAdditionFSM(vam, initializer, sampler, rng)

This puts together the model and a sampler. This is what we think of as
a simulation. We will organize this like it's a finite state machine,
so it will have an initializer, a dynamics, and an observer.
"""
mutable struct VectorAdditionFSM{
        VasT <: VectorAdditionSystem, 
        VasS <: VectorAdditionState,
        S,
        R <: AbstractRNG,
        I <: Function,
        O <: Function
    }  
    # This is rules about the simulation. It's part of the dynamics.
    vas::VasT
    state::VasS
    # The sampler does hold state of which events are enabled.
    sampler::S
    # The random number generator has state, too.
    rng::R

    # The initializer, aka iota.
    initializer::I
    is_initialized::Bool

    # Observer, aka lambda.
    observer::O
end


"""
This creates a simulation, taking as input a model, an initializer, a sampler,
and a random number generator. This combines several steps we could do
separately.

 1. Create the rules for the simulation, the VectorAdditionSystem.
 2. Combine the immutable rules with mutable state into a VectorAdditionModel.
 3. 
"""
function VectorAdditionFSM(vas, initializer, sampler, rng)
    state = VectorAdditionState(zero_state(vas), 0.0)
    VectorAdditionFSM(vas, state, sampler, rng, initializer, false, vam_event_observer)
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
        fire!(fsm.sampler, fsm.vas, fsm.state.state, fsm.initializer, fsm.state.when, fsm.rng)
        fsm.is_initialized = true
    end
    (when, what) = next(fsm.sampler, fsm.state.when, fsm.rng)
    if when < Inf
        action = vas_delta(fsm.vas, what)
        fsm.state.when = when
        fire!(fsm.sampler, fsm.vas, fsm.state.state, action, fsm.state.when, fsm.rng)
        fsm.observer(fsm.state, when, what)
    else
        (when, what)
    end
end

end # module VectorAddition
