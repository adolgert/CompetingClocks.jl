# # Sample Main Loop

# ## Introduction
# Let's walk through a small simulation so that we can see how CompetingClocks
# could appear in the main loop.
# This models individuals wandering around on a checkerboard, where no
# two individuals can occupy the same space. First, let's import libraries.

import Base
using Distributions
using Random
using SparseArrays
using Test
using CompetingClocks

# There is a physical state for the model, separate from the state of
# the bag of clocks and separate from the time. Most of the board is empty,
# so let's use a spare matrix to represent the locations of individuals.

mutable struct PhysicalState
    board::SparseMatrixCSC{Int64, Int64}
end

## They can move in any of four directions.
@enum Direction Up Left Down Right
const DirectionDelta = Dict(
    Up => CartesianIndex(-1, 0),
    Left => CartesianIndex(0, -1),
    Down => CartesianIndex(1, 0),
    Right => CartesianIndex(0, 1),
    );

# The simulation, itself, carries the state of the clocks in the sampler, as
# well as the physical state. We'll call it a finite state machine (FSM)
# because it has the traits of a
# [Moore machine](https://en.wikipedia.org/wiki/Moore_machine).

mutable struct SimulationFSM{Sampler}
    physical::PhysicalState
    sampler::Sampler
    when::Float64
    rng::Xoshiro
end

# ## Main Loop
# The main loop will ask what event happens next to the state. When that
# event changes the state, the loop will update the set of possible next
# events by disabling outdated ones and enabling new ones. The calls to
# CompetingClocks are:
#
# * `next(sampler, current time, random number generator (RNG))`
# * `enable!(sampler, event ID, distribution, current time, start time of clock, RNG)`
# * `disable!(sampler, event ID, current time)`
#
# There are a lot of samplers in CompetingClocks to choose from. This example uses `CombinedNextReaction`
# algorithm, which has good performance for a variety of distributions. Samplers in CompetingClocks
# require two type parameters, a key type for clocks and the type used to represent time.
# In this case, the clock key type fully represents an event, giving the ID of the individual,
# where they start, and which direction they may move.

const ClockKey = Tuple{Int,CartesianIndex{2},Direction}

function run(event_count)
    Sampler = CombinedNextReaction{ClockKey,Float64}
    physical = PhysicalState(zeros(Int, 10, 10))
    @test showable(MIME("text/plain"), physical)
    sim = SimulationFSM{Sampler}(
        physical,
        Sampler(),
        0.0,
        Xoshiro(2947223)
    )
    initialize!(sim.physical, 9, sim.rng)
    current_events = allowed_moves(sim.physical)
    for event_id in current_events
        enable!(sim.sampler, event_id, Weibull(1.0), 0.0, 0.0, sim.rng)
    end

    for i in 1:event_count
        (when, what) = next(sim.sampler, sim.when, sim.rng)
        if isfinite(when) && !isnothing(what)
            sim.when = when
            move!(sim.physical, what)
            next_events = allowed_moves(sim.physical)
            for remove_event in setdiff(current_events, next_events)
                disable!(sim.sampler, remove_event, when)
            end
            for add_event in setdiff(next_events, current_events)
                enable!(sim.sampler, add_event, Weibull(1.0), when, when, sim.rng)
            end
            current_events = next_events
            @show (when, what)
        end
    end
end;

# For this checkerboard with wandering individuals, the allowed moves are
# moves by any individual to any free square on the board. The set of allowed
# moves is precisely the set of enabled clocks, so it stores `ClockKey`s.

function allowed_moves(physical::PhysicalState)
    allowed = Set{ClockKey}()
    row, col, value = findnz(physical.board)
    for ind_idx in eachindex(value)
        location = CartesianIndex((row[ind_idx], col[ind_idx]))
        for (direction, offset) in DirectionDelta
            if checkbounds(Bool, physical.board, location + offset)
                if physical.board[location + offset] == 0
                    push!(allowed, (value[ind_idx], location, direction))
                end
            end
        end
    end
    return allowed
end;

# ## Changes to the state of the board
# We set up the board with an initializer that places individuals at random.
# We move one individual at a time, when their next event happens.

function initialize!(physical::PhysicalState, individuals::Int, rng)
    physical.board .= 0
    dropzeros!(physical.board)
    locations = zeros(CartesianIndex{2}, individuals)
    for ind_idx in 1:individuals
        loc = rand(rng, CartesianIndices(physical.board))
        while physical.board[loc] != 0
            loc = rand(rng, CartesianIndices(physical.board))
        end
        locations[ind_idx] = loc
        physical.board[loc] = ind_idx
    end
end;


function move!(physical::PhysicalState, event_id)
    (individual, previous_location, direction) = event_id
    next_location = previous_location + DirectionDelta[direction]
    ## This sets the previous board value to zero.
    SparseArrays.dropstored!(physical.board, previous_location.I...)
    physical.board[next_location] = individual
end;

# ## How it runs
# A run of this simulation produces a sequence of moves, no two happening
# at the same time because this is a continuous-time simulation.

run(10)
