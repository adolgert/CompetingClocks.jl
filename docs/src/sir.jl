# # Non-Markovian SIR Model

using Random
using Plots
using Distributions
using Fleck

# In this tutorial we demonstrate construction of a SIR (susceptible-infectious-removed) model
# with non-exponential recovery times. Infection events occur at the points of a Poisson process,
# which is equivalent to using a single exponential clock whose rate corresponds to the overall
# rate of infection in the population. The infection rate is ``\beta c S I / N`` where ``N`` is
# the total population size, making this a *frequency-dependent* force of infection term rather
# than pure mass action (although as ``N`` does not change in this example, the difference is moot).
# Clocks for recovery events follow an arbitrary distribution.

# ## Model structure

# We use a struct that stores the type of the recovery distribution as a type parameter.
# Additionally, because recovery clocks need unique keys, we define the method `get_key!`
# which retrives an integer key and increments the stored key counter.

# The `initialize!` method uses the initial state to enable the infection clock and
# recovery clocks for each initial infectious person. Note that the keys are a tuple
# where the first element is a `Symbol` giving the event type, and the second element
# is an integer.

mutable struct SIRNonMarkov{T<:Distribution}
    state::Vector{Int}
    parameters::Vector{Float64}
    next_key::Int
    recovery_distribution::T
    time::Float64
end

function get_key!(model::SIRNonMarkov)
    key = model.next_key
    model.next_key += 1
    return key
end

function initialize!(model::SIRNonMarkov, sampler, rng)
    (β, c, γ) = model.parameters
    ## enable the infection clock
    enable!(sampler, (:infection, get_key!(model)), Exponential(1.0/(β*c*model.state[2]/sum(model.state)*model.state[1])), model.time, model.time, rng)
    ## enable the recovery clocks
    for _ in 1:model.state[2]
        enable!(sampler, (:recovery, get_key!(model)), model.recovery_distribution, model.time, model.time, rng)
    end
end;

# ## Model update

# In the model update function, we use the first element of the clock key to determine the 
# event type corresponding to the clock that fired, and apply the corresponding logic.
# Infection events disable and enable the infection event with a new rate, and enable a
# recovery clock for the newly infectious individual. Recovery events simply disable the 
# clock associated to that event. Both events update the state vector.

function step!(model::SIRNonMarkov, sampler::SSA{K,T}, when::T, which::K, rng) where {K,T}
    (β, c, γ) = model.parameters
    model.time = when
    if first(which) == :infection
        model.state[1] -= 1
        model.state[2] += 1
        ## disable and reenable the infection clock after accounting for the new rate
        disable!(sampler, which, model.time)
        enable!(sampler, which, Exponential(1.0/(β*c*model.state[2]/sum(model.state)*model.state[1])), model.time, model.time, rng)
        ## enable a recovery event for the newly infected person
        enable!(sampler, (:recovery, get_key!(model)), model.recovery_distribution, model.time, model.time, rng)
    elseif first(which) == :recovery
        model.state[2] -= 1
        model.state[3] += 1
        disable!(sampler, which, model.time)
    else
        error("unrecognized clock key: $(which)")
    end
end;

# ## Simulation

# We first set parameters.

tmax = 40.0
initial_state = [990, 10, 0]
p = [0.05, 10.0, 4.0] 

seed = 456959517
rng = MersenneTwister(seed);

# Next we generate the model struct and the sampler object. Here
# we choose the `CombinedNextReaction` sampler type.

sirmodel = SIRNonMarkov(deepcopy(initial_state), p, 0, Dirac(p[3]), 0.0)
sampler = CombinedNextReaction{Tuple{Symbol,Int},Float64}();

# Now we may intiialize and run the model.
# We preallocate a matrix to store model output. Note that in the simple
# SIR model with only infection and recovery events, a maximum of ``S I + I``
# events is possible.

function run_sir!(model, sampler, tmax, rng)

    output = zeros(prod(sirmodel.state[1:2])+sum(sirmodel.state[1:2])+1, 4)
    nout = 1
    output[nout,:] = [sirmodel.time; sirmodel.state]
    nout += 1;

    initialize!(model, sampler, rng)

    (when, which) = next(sampler, model.time, rng)

    while when <= tmax
        step!(model, sampler, when, which, rng)
        (when, which) = next(sampler, model.time, rng)
    
        output[nout,:] .= [model.time; model.state]
        nout += 1
    end
    
    return output, nout
end

(output, nout) = run_sir!(sirmodel, sampler, tmax, rng);

# Finally we can plot the sampled trajectory.

plot(
    output[1:nout-1,1], 
    output[1:nout-1,2:end],
    label=["S" "I" "R"],
    xlabel="Time",
    ylabel="Number"
)
