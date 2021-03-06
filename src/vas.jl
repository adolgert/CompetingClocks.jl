# A Vector Addition System
import Distributions


"""
A `VectorAdditionSystem` is an older form of simulation that is a lot
like Petri nets. The state is a vector.
The system is a list of transitions. Each transition
is an array of values. Negative numbers mean the transition needs to
take this many tokens from the state, meaning the state at those indices
must be an integer at least that large. Positive numbers mean the
transition places tokens into the state. Unlike chemical simulations,
the rate doesn't depend on the number of combinations of species present.
"""
struct VectorAdditionSystem
    take::Array{Int, 2}  # states x transitions
    give::Array{Int, 2}  # states x transitions
    rates::Vector{Function}  # length is transitions
end


function zero_state(vas::VectorAdditionSystem)
    zeros(Int, size(vas.take, 1))
end


function vas_input(vas::VectorAdditionSystem, transition_idx)
    state_change = vas.give[:, transition_idx] - vas.take[:, transition_idx]
    let delta = state_change
        state -> begin state .+= state_change end
    end
end


function vas_initial(vas::VectorAdditionSystem, initial_state)
    let initial = initial_state
        state -> begin state .= initial end
    end
end


function fire!(visitor::Function, vas::VectorAdditionSystem, state, modify_state, rng)
    state_prime = copy(state)
    modify_state(state_prime)
    for rate_idx in eachindex(vas.rates)
        was_enabled = all(state .- vas.take[:, rate_idx] .>= 0)
        now_enabled =  all(state_prime .- vas.take[:, rate_idx] .>= 0)
        if was_enabled && !now_enabled
            visitor(rate_idx, Distributions.Exponential(1), :Disabled, rng)
        elseif !was_enabled && now_enabled
            visitor(rate_idx, vas.rates[rate_idx](state_prime), :Enabled, rng)
        elseif was_enabled && now_enabled
            ratefunc = vas.rates[rate_idx]
            former_rate = ratefunc(state)
            current_rate = ratefunc(state_prime)
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
    sampler::DirectCall{Int}
end


function simstep!(fsm::VectorAdditionFSM, state_update::Function, rng::AbstractRNG)
    visitor = (clock, dist, enabled, gen) -> begin
        set_clock!(fsm.sampler, clock, dist, enabled, gen)
    end
    fire!(visitor, fsm.vam.vas, fsm.vam.state, state_update, rng)
    next(fsm.sampler, fsm.vam.when, rng)
end
