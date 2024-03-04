using DataStructures

export FirstToFire

"""
    FirstToFire{KeyType,TimeType}()

Construct a FirstToFire.
As soon as a distribution is enabled, this draws a value from the distribution.
The soonest to fire wins. When a clock is disabled, its future firing time is
removed from the list. There is no memory of previous firing times.
"""
struct FirstToFire{K,T} <: SSA{K,T}
    firing_queue::MutableBinaryMinHeap{OrderedSample{K,T}}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{K,Int}
end


function FirstToFire{K,T}() where {K,T}
    heap = MutableBinaryMinHeap{OrderedSample{K,T}}()
    state = Dict{K,Int}()
    FirstToFire{K,T}(heap, state)
end


# Finds the next one without removing it from the queue.
function next(propagator::FirstToFire{K,T}, when::T, rng::AbstractRNG) where {K,T}
    least = if !isempty(propagator.firing_queue)
        top(propagator.firing_queue)
    else
        OrderedSample(nothing, typemax(T))
    end
    @debug("FirstToFire.next queue length ",
            length(propagator.firing_queue), " least ", least)
    (least.time, least.key)
end


function enable!(
    propagator::FirstToFire{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T, rng::AbstractRNG) where {K,T}

    if te < when
        when_fire = te + rand(rng, truncated(distribution, when - te, typemax(T)))
    else
        when_fire = te + rand(rng, distribution)
    end
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock]
        update!(propagator.firing_queue, heap_handle, OrderedSample{K,T}(clock, when_fire))
    else
        heap_handle = push!(propagator.firing_queue, OrderedSample{K,T}(clock, when_fire))
        propagator.transition_entry[clock] = heap_handle
    end
end


function disable!(propagator::FirstToFire{K,T}, clock::K, when::T) where {K,T}
    heap_handle = propagator.transition_entry[clock]
    delete!(propagator.firing_queue, heap_handle)
    delete!(propagator.transition_entry, clock)
end
