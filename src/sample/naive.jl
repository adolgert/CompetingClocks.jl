using DataStructures
import Base: ==, <, >

export FirstToFire

"""
A record of a transition and the time.
It's sortable by time. Immutable.
"""
struct NRTransition{T}
	key::T
	time::Float64
end


function Base.:<(a::NRTransition, b::NRTransition)
	a.time < b.time
end


function Base.isless(a::NRTransition, b::NRTransition)
    return isless(a.time, b.time)
end


function Base.:>(a::NRTransition, b::NRTransition)
    a.time > b.time
end


function ==(a::NRTransition, b::NRTransition)
    a.time == b.time
end


"""
Every time a clock is enabled, this draws a putative time at which it would
fire. The next to fire is chosen, and those that are disabled are removed from
the list.
"""
struct FirstToFire{T}
    firing_queue::MutableBinaryHeap{NRTransition{T}}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{T,Int}
    disabled::Set{T}
end


"""
Construct a FirstToFire.
This doesn't require a clock that remembers integrated hazard.
This sampler is inappropriate if any transition is either
re-enabled after being disabled or has its distribution
modified while enabled. It's OK if each transition fires
only once.

BUT I can't create a model where this sampler's output
varies from First Reaction. If anyone can show me when
this fails to work, I'd be grateful. Even the sis.jl example works.
"""
function FirstToFire(KeyType::Type)
    heap = MutableBinaryMinHeap{NRTransition{KeyType}}()
    state = Dict{KeyType,Int}()
    FirstToFire{KeyType}(heap, state, Set{KeyType}())
end


# Finds the next one without removing it from the queue.
function next(propagator::FirstToFire, when::Float64, rng::AbstractRNG)
    least = if !isempty(propagator.firing_queue)
        top(propagator.firing_queue)
    else
        NRTransition(nothing, Inf)
    end
    @debug("FirstToFire.next queue length ",
            length(propagator.firing_queue), " least ", least)
    (least.time, least.key)
end


function enable!(
    propagator::FirstToFire{T}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    when_fire = rand(rng, truncated(distribution, when - te, Inf))
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock]
        update!(propagator.firing_queue, heap_handle, NRTransition{T}(clock, when_fire))
    else
        heap_handle = push!(propagator.firing_queue, NRTransition{T}(clock, when_fire))
        propagator.transition_entry[clock] = heap_handle
    end
end


function disable!(propagator::FirstToFire{T}, clock::T, when::Float64) where {T}
    heap_handle=propagator.transition_entry[clock]
    delete!(propagator.firing_queue, heap_handle)
    delete!(propagator.transition_entry, clock)
end
