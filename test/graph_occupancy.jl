using Distributions
using CompetingClocks
using Graphs
using Random
using StatsBase
using Base

const Time = Float64

"""
A transition describes how a particle can move from one node to another along
an edge. The distribution is the distribution of times for the jump, as a competing
process. The relative enabling time modifies that distribution to start a little
before or after the current time, if we wish. If the transition has memory, then
when a particle arrives at the node, it keeps the previous enabling time and re-enables
the transition. If it has memory, it will store the duration of time it has already
been enabled.
"""
mutable struct Transition
    distribution::UnivariateDistribution
    relative_enabling_time::Time
    memory::Bool
    consumed::Time
end


"""
This is the struct that describes both the model rules (graph shape and
exact transitions on that graph) and the state of the system (location of the
particle and current time).
"""
mutable struct GraphOccupancy
    g::SimpleGraph{Int64}  # Graph of nodes and edges to other nodes.
    transition::Dict{Tuple{Int,Int},Transition}
    vertex::Int64  # Current node on which particle resides.
    when::Time
end


"""
This is the type of the key the model uses to identify transitions.
Here, it is a tuple of (node from which to hop, node to which to hop).
"""
keyspace(::Type{GraphOccupancy}) = Tuple{Int,Int}
Base.length(go::GraphOccupancy) = nv(go.g)


"""
Make a model. The rng is a random number generator. Initialize that generator
with a different seed to get a different system. The features is a dictionary
of boolean values to tell this model how complicated it should be. We want
to turn features on/off in order to ascertain when and why different samplers
disagree about how to simulate the model.
"""
function GraphOccupancy(features::AbstractDict{String,Bool}, rng::AbstractRNG)
    # Make a random graph with random transitions.
    seed = rand(rng, 1:100000)
    g = barabasi_albert(8, 4, 3, seed=seed)  # add 4 nodes, each connected to 3 existing.
    # It is some work to make random transitions. An actual model would be simpler.
    transition = Dict{Tuple{Int,Int},Transition}()
    for e in edges(g)
        for (l, r) in [(src(e), dst(e)), (dst(e), src(e))]
            distributions = [
                Exponential(0.8 + 0.4 * rand(rng)),
            ]
            if features["weibull"]
                push!(distributions, Weibull(rand(rng, [1, 1.5, 2])))
            end
            if features["gamma"]
                push!(distributions, Gamma(
                    rand(rng, [1, 2, 3, 5]),
                    rand(rng, [1, 2])
                ))
            end
            distribution = sample(rng, distributions)
            delta = 0.1 + rand(rng)
            reltimes = [0]  # This is a shift to the start of the distribution.
            relweights = [0.7]  # The probability of choosing this shift value.
            if features["past"]
                push!(reltimes, -delta)
                push!(relweights, 0.2)
            end
            if features["future"]
                push!(reltimes, delta)
                push!(relweights, 0.1)
            end
            relative_enabling_time = sample(rng, reltimes, Weights(relweights))
            if features["memory"]
                memory = sample(rng, [false, true], Weights([0.7, 0.3]))
            else
                memory = false
            end
            transition[(l, r)] = Transition(distribution, relative_enabling_time, memory, 0.0)
        end
    end
    vertex = 1
    when = 0.0
    GraphOccupancy(g, transition, vertex, when)
end


"""
Before the simulation starts, the initial state should mean some transitions are
enabled.
"""
function initial_enabling(go::GraphOccupancy, sampler, rng::AbstractRNG)
    vertex = go.vertex
    for neighbor in neighbors(go.g, vertex)
        trans = go.transition[(vertex, neighbor)]
        enable!(sampler, (vertex, neighbor), trans.distribution,
            go.when + trans.relative_enabling_time, go.when, rng)
    end
end


function step!(go::GraphOccupancy, sampler, when::Time,
        which::Tuple{Int64,Int64}, rng::AbstractRNG
        )
    if go.vertex != which[1]
        println("Expected vertex to match $(go.vertex) but got $which")
    end
    if go.when >= when
        println("Expected old time $(go.when) to be less than new time $when")
        debug_trans = go.transition[which]
        println("Transition was $(debug_trans.distribution)")
        println("Transition shift $(debug_trans.relative_enabling_time)")
        println("Transition memory $(debug_trans.memory)")
        println("Transition consumed $(debug_trans.consumed)")
    end
    @assert go.vertex == which[1]
    @assert which[2] ∈ neighbors(go.g, which[1])
    @assert go.when <= when

    # First disable previously-enabled hops from one node to the next.
    for dis_neighbor in neighbors(go.g, go.vertex)
        trans = go.transition[(go.vertex, dis_neighbor)]
        if trans.memory
            if (go.vertex, dis_neighbor) == which
                # The transition that fires has its memory reset.
                trans.consumed = 0.0
            else
                # A transition with memory will increase its hazard rate by
                # time-shifting the next time it's enabled.
                trans.consumed += when - go.when
            end
        end
        disable!(sampler, (go.vertex, dis_neighbor), when)
    end

    resident = go.vertex
    duration = when - go.when

    go.vertex = which[2]  # The transition key is (from vertex, to vertex).
    go.when = when

    # Then enable the new ones that start at the next node.
    for en_neighbor in neighbors(go.g, go.vertex)
        trans = go.transition[(go.vertex, en_neighbor)]
        if trans.memory
            # The memory means this transition says its enabling time is in the past.
            enable!(sampler, (go.vertex, en_neighbor), trans.distribution,
                go.when + trans.relative_enabling_time - trans.consumed, go.when, rng)
        else
            enable!(sampler, (go.vertex, en_neighbor), trans.distribution,
                go.when + trans.relative_enabling_time, go.when, rng)
        end
    end

    # Think of this as an observer of the state of the model.
    (resident, duration)
end


"""
This records how often a transition fired and for how long.
It's used for watching output of the simulation.
"""
mutable struct TransitionObserver
    min_duration::Time
    max_duration::Time
    mean_duration::Time
    total_duration::Time
    call_cnt::Int64
end


"""
This is the test. Run this with a sampler and check the occupancy numbers.
"""
function run_graph_occupancy(groc::GraphOccupancy, finish_time, sampler, rng)
    initial_enabling(groc, sampler, rng)
    occupancy = zeros(Time, length(groc))
    KeyType = keyspace(GraphOccupancy)
    observations = Dict{KeyType,TransitionObserver}()

    when, which = next(sampler, groc.when, rng)
    while which !== nothing && when < finish_time
        resident_node, duration = step!(groc, sampler, when, which, rng)

        # The stopping time is finish_time, so don't count beyond the stopping time.
        if when >= finish_time
            duration -= (when - finish_time)
            @assert duration > 0.0
        end
        occupancy[resident_node] += duration

        if !haskey(observations, which)
            observations[which] = TransitionObserver(Inf, 0.0, 0.0, 0.0, 0)
        end
        observations[which].total_duration += duration
        observations[which].min_duration = min(duration, observations[which].min_duration)
        observations[which].max_duration = max(duration, observations[which].max_duration)
        observations[which].call_cnt += 1

        when, which = next(sampler, groc.when, rng)
    end
    for transobv in values(observations)
        transobv.mean_duration = transobv.total_duration / transobv.call_cnt
    end
    occupancy, observations
end


function sample_run_graph_occupancy(options=Set{String}())
    model_rng = Xoshiro(349827)
    features = Dict{String,Bool}(
        "weibull"=>false, "gamma"=>false, "past"=>false, "future"=>false,
        "memory"=>false)
    for turn_on in options
        if turn_on ∉ keys(features)
            println("Trying to turn on $turn_on but it isn't a feature in $features")
        end
        features[turn_on] = true
    end
    single_groc = GraphOccupancy(features, model_rng)

    groc = deepcopy(single_groc)
    sampler = ChatReaction{keyspace(GraphOccupancy)}()
    rng = Xoshiro(2342374)
    occupancy = run_graph_occupancy(groc, 1e6, sampler, rng)
    return occupancy
end
