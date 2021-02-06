using DataStructures

"""
This implements an incorrect algorithm which is the most common
intuitive choice. Each trajectory's distribution looks correct,
but the master equation isn't correct.
"""

struct NaiveSampler
    firing_queue::MutableBinaryHeap{NRTransition,Base.Order.ForwardOrdering}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{Any,Int}
    disabled::Set{Any}
    init::Bool
end


"""
Construct a NaiveSampler.
This doesn't require a clock that remembers integrated hazard.
This sampler is inappropriate if any transition is either
re-enabled after being disabled or has its distribution
modified while enabled. It's OK if each transition fires
only once.

BUT I can't create a model where this sampler's output
varies from First Reaction. If anyone can show me when
this fails to work, I'd be grateful. Even the sis.jl example works.
"""
function NaiveSampler()
    heap=mutable_binary_minheap(NRTransition)
    state=Dict{Any,Int}()
    NaiveSampler(heap, state, Set{Any}(), true)
end


# Finds the next one without removing it from the queue.
function Next(propagator::NaiveSampler, system, rng)
    if propagator.init
        Hazards(system, rng) do clock, now, updated, rng2
            NaiveObserve(propagator, clock, now, updated, rng2)
        end
        propagator.init=false
    end

    NotFound=NRTransition(nothing, Inf)
    if !isempty(propagator.firing_queue)
        least=top(propagator.firing_queue)
    else
        least=NotFound
    end
    @debug("SampleSemiMarkov.next queue length ",
            length(propagator.firing_queue), " least ", least)
    (least.time, least.key)
end


"""
Returns an observer of intensities to decide what to
do when they change.
"""
function Observer(propagator::NaiveSampler)
    function nrobserve(clock, time, updated, rng)
        NaiveObserve(propagator, clock, time, updated, rng)
    end
end


function NaiveObserve(propagator::NaiveSampler, clock,
            time, updated, rng)
    key=clock
    if updated==:Fired || updated==:Disabled
        heap_handle=propagator.transition_entry[key]
        # We store distributions in order to calculate remaining hazard
        # which will happen AFTER the state has changed.
        update!(propagator.firing_queue, heap_handle,
            NRTransition(key, -1.))
        todelete=pop!(propagator.firing_queue)
        delete!(propagator.transition_entry, key)
        push!(propagator.disabled, clock)

    elseif updated==:Enabled
        # if haskey(propagator.disabled, clock)
        #     error("Cannot re-enable a transition with this sampler.")
        # end
        when_fire=Sample(clock.intensity, time, rng)
        heap_handle=push!(propagator.firing_queue,
                NRTransition(key, when_fire))
        propagator.transition_entry[key]=heap_handle

    elseif updated==:Modified
        # if haskey(propagator.disabled, clock)
        #     error("Cannot modify a transition with this sampler.")
        # end
        when_fire=Sample(clock.intensity, time, rng)
        heap_handle=propagator.transition_entry[key]
        update!(propagator.firing_queue, heap_handle,
                NRTransition(key, when_fire))
    else
        assert(false)
    end
end
