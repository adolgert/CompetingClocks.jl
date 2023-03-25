using DataStructures

export FirstToFire

"""
    FirstToFire(KeyType::Type)

Construct a FirstToFire.
As soon as a distribution is enabled, this draws a value from the distribution.
The soonest to fire wins. When a clock is disabled, its future firing time is
removed from the list. There is no memory of previous firing times.
"""
struct FirstToFire{T}
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{T,Int}
end


function FirstToFire(KeyType::Type)
    heap = MutableBinaryMinHeap{OrderedSample{KeyType}}()
    state = Dict{KeyType,Int}()
    FirstToFire{KeyType}(heap, state)
end


# Finds the next one without removing it from the queue.
function next(propagator::FirstToFire, when::Float64, rng::AbstractRNG)
    least = if !isempty(propagator.firing_queue)
        top(propagator.firing_queue)
    else
        OrderedSample(nothing, Inf)
    end
    @debug("FirstToFire.next queue length ",
            length(propagator.firing_queue), " least ", least)
    (least.time, least.key)
end


function enable!(
    propagator::FirstToFire{T}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    if te < when
        when_fire = te + rand(rng, truncated(distribution, when - te, Inf))
    else
        when_fire = te + rand(rng, distribution)
    end
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock]
        update!(propagator.firing_queue, heap_handle, OrderedSample{T}(clock, when_fire))
    else
        heap_handle = push!(propagator.firing_queue, OrderedSample{T}(clock, when_fire))
        propagator.transition_entry[clock] = heap_handle
    end
end


function disable!(propagator::FirstToFire{T}, clock::T, when::Float64) where {T}
    heap_handle = propagator.transition_entry[clock]
    delete!(propagator.firing_queue, heap_handle)
    delete!(propagator.transition_entry, clock)
end
