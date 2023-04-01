using Distributions
using Fleck
using Graphs
using Random
using StatsBase
using Base

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
    relative_enabling_time::Float64
    memory::Bool
    consumed::Float64
end


mutable struct GraphOccupancy
    g::SimpleGraph{Int64}
    transition::Dict{Tuple{Int,Int},Transition}
    vertex::Int64
    when::Float64
end


"""
This is the type of the key the model uses to identify transitions.
Here, it is a (Int, Int) tuple for GraphOccupancy.
"""
keyspace(::Type{GraphOccupancy}) = Tuple{Int,Int}
Base.length(go::GraphOccupancy) = nv(go.g)


function GraphOccupancy(rng::AbstractRNG)
    # Make a random graph.
    g = barabasi_albert(8, 4, 3, seed=123)  # add 4 nodes, each connected to 3 existing.
    # It is some work to make random transitions. An actual model would be simpler.
    transition = Dict{Tuple{Int,Int},Transition}()
    for e in edges(g)
        for (l, r) in [(src(e), dst(e)), (dst(e), src(e))]
            distributions = [
                Exponential(0.8 + 0.4 * rand(rng)),
                Weibull(rand(rng, [0.5, 1, 1.5, 2])),
                Gamma(
                    rand(rng, [1, 2, 3, 5]),
                    rand(rng, [1, 2])
                )
            ]
            distribution = sample(rng, distributions)
            delta = rand(rng)
            relative_enabling_time = sample(rng, [-delta, 0, delta], Weights([0.2, 0.7, 0.1]))
            memory = sample(rng, [false, true], Weights([0.7, 0.3]))
            transition[(l, r)] = Transition(distribution, relative_enabling_time, memory, 0.0)
        end
    end
    vertex = 1
    when = -5.0 + 10 * rand(rng)  # Why not start with a weird time?
    GraphOccupancy(g, transition, 1, when)
end


function initial_enabling(go::GraphOccupancy, sampler, rng::AbstractRNG)
    vertex = go.vertex
    for neighbor in neighbors(go.g, vertex)
        trans = go.transition[(vertex, neighbor)]
        enable!(sampler, (vertex, neighbor), trans.distribution,
            go.when + trans.relative_enabling_time, go.when, rng)
    end
end


function step!(go::GraphOccupancy, sampler, when::Float64,
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
    @assert go.when < when

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
This is the test. Run this with a sampler and check the occupancy numbers.
"""
function run_graph_occupancy(groc::GraphOccupancy, step_cnt, sampler, rng)
    initial_enabling(groc, sampler, rng)
    occupancy = zeros(Float64, length(groc))
    for step_idx in 1:step_cnt
        when, which = next(sampler, groc.when, rng)
        resident_node, duration = step!(groc, sampler, when, which, rng)
        occupancy[resident_node] += duration
    end
    occupancy
end


function sample_run_graph_occupancy()
    model_rng = Xoshiro(349827)
    single_groc = GraphOccupancy(model_rng)

    groc = deepcopy(single_groc)
    sampler = ChatReaction{keyspace(GraphOccupancy)}()
    rng = Xoshiro(2342374)
    occupancy = run_graph_occupancy(groc, 1000, sampler, rng)
    return occupancy
end
