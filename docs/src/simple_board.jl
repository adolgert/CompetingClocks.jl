# # Sample Main Loop

# ## Introduction
# Let's walk through a small simulation so that we can see how CompetingClocks
# could appear in the main loop.
# This models individuals wandering around on a checkerboard, where no
# two individuals can occupy the same space. First, let's import libraries.

import Base
using Distributions
using Random
using Test
using CompetingClocks

# There is a physical state for the model, separate from the state of
# the bag of clocks and separate from the time. Most of the board is empty,
# so let's use a spare matrix to represent the locations of individuals.

mutable struct PhysicalState
    board::Matrix{Int64}
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
# well as the physical state. We'll call it a
# [finite state machine](https://en.wikipedia.org/wiki/Moore_machine).

mutable struct SimulationFSM{Sampler}
    physical::PhysicalState
    sampler::Sampler
end

# ## Main Loop
# The main loop will ask what event happens next to the state. When that
# event changes the state, the loop will update the set of possible next
# events by disabling outdated ones and enabling new ones. The calls to
# CompetingClocks are:
#
# * `next(sampler)`
# * `enable!(sampler, event ID, distribution)`
# * `disable!(sampler, event ID)`
#
# There are a lot of samplers in CompetingClocks to choose from. Samplers in CompetingClocks
# require two type parameters, a key type for clocks and the type used to represent time.
# In this case, the clock key type fully represents an event, giving the ID of the individual,
# where they start, and which direction they may move.

struct ClockKey
    individual::Int
    location::CartesianIndex{2}
    direction::Direction
end

function run(event_count)
    rng = Xoshiro(2947223)
    builder = SamplerBuilder(ClockKey, Float64)
    sampler = SamplingContext(builder, rng)
    physical = PhysicalState(zeros(Int, 10, 10))
    sim = SimulationFSM(
        physical,
        sampler,
    )
    initialize!(sim.physical, 9, rng)
    current_events = allowed_moves(sim.physical)
    for event_id in current_events
        enable!(sim.sampler, event_id, Weibull(1.0))
    end

    for i in 1:event_count
        (when, what) = next(sim.sampler)
        if isfinite(when) && !isnothing(what)
            fire!(sim.sampler, what, when)
            delete!(current_events, what)
            move!(sim.physical, what)
            next_events = allowed_moves(sim.physical)
            update_events!(sampler, current_events, next_events, Weibull(1.0))
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
    occupied_locations = findall(!iszero, physical.board)
    for location in occupied_locations
        individual = physical.board[location]
        for (direction, offset) in DirectionDelta
            if checkbounds(Bool, physical.board, location + offset)
                if physical.board[location + offset] == 0
                    push!(allowed, ClockKey(individual, location, direction))
                end
            end
        end
    end
    return allowed
end;

# In practice, you would track dependencies among individuals in order to
# accelerate a simulation rather than checking all events in `allowed_moves`. 
# The next function updates the current events in the sampler.

function update_events!(sampler, old_events, new_events, distribution)
    for remove in setdiff(old_events, new_events)
        disable!(sampler, remove)
    end
    for add in setdiff(new_events, old_events)
        enable!(sampler, add, distribution)
    end
end;

# ## Changes to the state of the board
# We set up the board with an initializer that places individuals at random.
# We move one individual at a time, when their next event happens.

function initialize!(physical::PhysicalState, individuals::Int, rng)
    physical.board .= 0
    for ind_idx in 1:individuals
        loc = rand(rng, CartesianIndices(physical.board))
        while physical.board[loc] != 0
            loc = rand(rng, CartesianIndices(physical.board))
        end
        physical.board[loc] = ind_idx
    end
end;


function move!(physical::PhysicalState, event::ClockKey)
    next_location = event.location + DirectionDelta[event.direction]
    physical.board[event.location] = 0
    physical.board[next_location] = event.individual
end;

# ## How it runs
# A run of this simulation produces a sequence of moves, no two happening
# at the same time because this is a continuous-time simulation.

run(10)
