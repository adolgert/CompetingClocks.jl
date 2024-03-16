# # Sample Main Loop

# Let's walk through a small simulation so that we can see how and when Fleck
# could appear in the main loop. We're going to split this up here in a
# way that is more simplistic than the Fleck interface requires.

import Base
using Distributions
using Random
using Test
using Fleck

# There is a physical state for the model, separate from the state of
# the bag of clocks and separate from the time.

mutable struct PhysicalState
    board::Matrix{Int}
end

function Base.show(io::IO, physical::PhysicalState)
    for row_idx in 1:size(physical.board, 1)
        println(io, join([string(x) for x in physical.board[row_idx, :]], ""))
    end
end


@enum Direction Up Left Down Right
const DirectionDelta = Dict(
    Up => CartesianIndex(-1, 0),
    Left => CartesianIndex(0, -1),
    Down => CartesianIndex(1, 0),
    Right => CartesianIndex(0, 1),
    )
const OppositeDirection = Dict(
    Up => Down,
    Down => Up,
    Left => Right,
    Right => Left
)

function over_neighbors(f::Function, physical::PhysicalState, location)
    for (direction, offset) in DirectionDelta
        if checkbounds(Bool, physical.board, location + offset)
            f(direction, location + offset)
        end
    end
end


## Use separate initializer in case you want to re-initialize the same memory.
function initialize!(physical::PhysicalState, sampler, individuals::Int, rng)
    physical.board .= 0
    locations = zeros(CartesianIndex{2}, individuals)
    for ind_idx in 1:individuals
        loc = rand(rng, CartesianIndices(physical.board))
        while physical.board[loc] != 0
            loc = rand(rng, CartesianIndices(physical.board))
        end
        locations[ind_idx] = loc
        physical.board[loc] = ind_idx
    end
    # During initialization, tell the sampler each individual can move to
    # any blank spot next to them.
    for mover_idx in 1:individuals
        location = locations[mover_idx]
        over_neighbors(physical, location) do direction, neighbor_location
            if physical.board[neighbor_location] == 0
                rate = Weibull(1.0)
                enable!(sampler, (mover_idx, direction), rate, 0.0, 0.0, rng)
            end
        end
    end
end


function move!(physical::PhysicalState, sampler, event_id, when, rng)
    (individual, direction) = event_id
    previous_location = findfirst(x -> x == individual, physical.board)
    next_location = previous_location + DirectionDelta[direction]

    # The individual who moves will disable all of its previous jumps.
    over_neighbors(physical, previous_location) do direction, neighbor_location
        if physical.board[neighbor_location] == 0
            disable!(sampler, (individual, direction), when)
        end
    end
    # If that individual had neighbors, those neighbors are now free to move
    # to the newly-vacated spot.
    over_neighbors(physical, previous_location) do direction, neighbor_location
        previous_neighbor = physical.board[neighbor_location]
        if previous_neighbor != 0
            rate = Weibull(1.0)
            now_move = OppositeDirection[direction]
            enable!(sampler, (previous_neighbor, now_move), rate, when, when, rng)
        end
    end

    physical.board[previous_location] = 0
    physical.board[next_location] = individual

    # The individual moved to a new spot that may block neighbor movements.
    over_neighbors(physical, next_location) do direction, neighbor_location
        new_neighbor = physical.board[neighbor_location]
        if new_neighbor != 0
            no_longer_move = OppositeDirection[direction]
            disable!(sampler, (new_neighbor, no_longer_move), when)
        end
    end
    # The individual now has a new set of jumps it can do, depending on
    # free spots nearby.
    over_neighbors(physical, next_location) do direction, neighbor_location
        if physical.board[neighbor_location] == 0
            rate = Weibull(1.0)
            enable!(sampler, (individual, direction), rate, when, when, rng)
        end
    end
end


# The simulation, itself, carries the state of the clocks, as well as the
# physical state. We'll call it a finite state machine (FSM) because it has
# the traits of a [Moore machine](https://en.wikipedia.org/wiki/Moore_machine).

mutable struct SimulationFSM{Sampler}
    physical::PhysicalState
    sampler::Sampler
    when::Float64
    rng::Xoshiro
end


function run()
    EventKey = Tuple{Int,Direction}
    TimeType = Float64
    Sampler = CombinedNextReaction{EventKey,TimeType}
    physical = PhysicalState(zeros(Int, 10, 10))
    @test showable(MIME("text/plain"), physical)
    sim = SimulationFSM{Sampler}(
        physical,
        Sampler(),
        0.0,
        Xoshiro(2947223)
    )
    initialize!(sim.physical, sim.sampler, 9, sim.rng)
    for i in 1:100
        (when, what) = next(sim.sampler, sim.when, sim.rng)
        if isfinite(when)
            sim.when = when
            move!(sim.physical, sim.sampler, what, sim.when, sim.rng)
            @show (when, what)
        else
            println("no possible moves")
            break
        end
    end
end

run()
