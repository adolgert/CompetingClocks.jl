# A Vector Addition System
import Distributions

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
"""
struct VectorAdditionSystem{T <: Function}
    take::Array{Int, 2}  # states x transitions
    give::Array{Int, 2}  # states x transitions
    rates::Vector{T}  # length is transitions
end


"""
    zero_state(vas::VectorAdditionSystem)

Return a state vector of all zeros.
"""
function zero_state(vas::VectorAdditionSystem)
    zeros(Int, size(vas.take, 1))
end


"""
    vas_delta(vas::VectorAdditionSystem, transition_idx)

Return a function taking a single argument `state` that applies the
state change assocaited with the transition indexed by `transition_idx`.
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
    fire!(visitor::Function, vas::VectorAdditionSystem, state, modify_state, rng)

Fire a transition. `visitor` is a function taking four arguments,  
`modify_state` is a function returned from `vas_delta`, or any other
function that accepts the state vector (argument `state`) as input and applies the update.
"""
function fire!(visitor::Function, vas::VectorAdditionSystem, state, modify_state, rng)
    former = copy(state)
    modify_state(state)
    for rate_idx in eachindex(vas.rates)
        was_enabled = all(former .- vas.take[:, rate_idx] .>= 0)
        now_enabled =  all(state .- vas.take[:, rate_idx] .>= 0)
        if was_enabled && !now_enabled
            visitor(rate_idx, Distributions.Exponential(1), :Disabled, rng)
        elseif !was_enabled && now_enabled
            visitor(rate_idx, vas.rates[rate_idx](state), :Enabled, rng)
        elseif was_enabled && now_enabled
            ratefunc = vas.rates[rate_idx]
            former_rate = ratefunc(former)
            current_rate = ratefunc(state)
            if former_rate != current_rate
                visitor(rate_idx, current_rate, :Changed, rng)
            end  # Else don't notify because rate is the same.
        end
    end
end


struct VectorAdditionModel
    vas::VectorAdditionSystem
    state::Vector{Int}
    when::Float64
end


using Random: AbstractRNG


struct VectorAdditionFSM
    vam::VectorAdditionModel
    sampler::Any
end


function simstep!(fsm::VectorAdditionFSM, state_update::Function, rng::AbstractRNG)
    visitor = (clock, dist, enabled, gen) -> begin
        set_clock!(fsm.sampler, clock, dist, enabled, gen)
    end
    fire!(visitor, fsm.vam.vas, fsm.vam.state, state_update, rng)
    next(fsm.sampler, fsm.vam.when, rng)
end
