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
    transitions::Array{Int, 2}  # states x transitions
    rates::Array{Float64, 1}  # length is transitions
end


function zero_state(vas::VectorAdditionSystem)
    zeros(Int, size(vas.transitions, 1))
end


function vas_input(vas::VectorAdditionSystem, transition_idx)
    state_change = vas.transitions[:, transition_idx]
    let delta = state_change
        state -> begin state .+= state_change end
    end
end


function vas_initial(vas::VectorAdditionSystem, initial_state)
    let initial = initial_state
        state -> begin state .= initial end
    end
end


function hazards!(visitor::Function, vas::VectorAdditionSystem, state, modify_state, rng)
    state_prime = copy(state)
    modify_state(state_prime)
    for rate_idx in eachindex(vas.rates)
        summed = vas.transitions[:, rate_idx] .+ state
        summed_prime = vas.transitions[:, rate_idx] .+ state_prime
        was_enabled = all(summed .>= 0)
        now_enabled = all(summed_prime .>= 0)
        if was_enabled && !now_enabled
            visitor(rate_idx, Distributions.Exponential(vas.rates[rate_idx]), :Disabled, rng)
        elseif !was_enabled && now_enabled
            visitor(rate_idx, Distributions.Exponential(vas.rates[rate_idx]), :Enabled, rng)
        end
    end
end
