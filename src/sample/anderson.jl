import DataStructures: MutableBinaryHeap, MutableBinaryMinHeap
using Random

"""
Next reaction by Hazards
Also called Anderson's method.
"""
struct TransitionRecord
	exponential_interval::Float64
	heap_handle::Int64
end


mutable struct NextReactionHazards
	firing_queue::MutableBinaryHeap{OrderedSample,Base.Order.ForwardOrdering}
	transition_state::Dict{Any,TransitionRecord}
	init::Bool
end


"""
Construct a Next Reaction sampler.
"""
function NextReactionHazards()
    heap = MutableBinaryMinHeap{OrderedSample}()
    @debug("SampleSemiMarkov.NextReactionHazards type ", typeof(heap))
    state = Dict{Any,TransitionRecord}()
    NextReactionHazards(heap, state, true)
end


# Finds the next one without removing it from the queue.
function Next(propagator::NextReactionHazards, system, rng)
	if propagator.init
		Hazards(system, rng) do clock, now, updated, rng2
			Enable(propagator, clock, now, updated, rng2)
	    end
	    propagator.init = false
	end

	NotFound = OrderedSample(nothing, Inf)
	if !isempty(propagator.firing_queue)
		least = top(propagator.firing_queue)
	else
		least = NotFound
	end
	@debug("SampleSemiMarkov.next queue length ",
			length(propagator.firing_queue), " least ", least)
	(least.time, least.key)
end


"""
Returns an observer of intensities to decide what to
do when they change.
"""
function Observer(propagator::NextReactionHazards)
	function nrobserve(clock, time, updated, rng)
		if updated == :Disabled || updated == :Fired
			Disable(propagator, clock, time, updated, rng)
		else
			Enable(propagator, clock, time, updated, rng)
		end
	end
end


function unit_hazard_interval(rng::MersenneTwister)
	-log(rand(rng))
end


# Enable or modify a hazard.
function Enable(propagator::NextReactionHazards, clock,
		now, updated, rng)
	key = clock
	clock_started = haskey(propagator.transition_state, key)
	if clock_started
		record = propagator.transition_state[key]
		when_fire = Putative(clock.intensity, now, record.exponential_interval)

		@assert(when_fire >= now)
		if record.heap_handle >= 0
			@debug("SampleSemiMarkov.enable keyu ", key, " interval ",
				record.exponential_interval, " when ", when_fire,
				" dist ", clock)
			update!(propagator.firing_queue, record.heap_handle,
				OrderedSample(key, when_fire))
		else
			record.heap_handle = push!(propagator.firing_queue,
				OrderedSample(key, when_fire))
			@debug("SampleSemiMarkov.enable keyp ", key, " interval ",
				record.exponential_interval, " when ", when_fire,
				" dist ", clock)
		end
	else
		firing_time, interval = MeasuredSample(clock.intensity, now, rng)
		@assert(firing_time >= now)
        handle = push!(propagator.firing_queue, OrderedSample(key, firing_time))
        @debug("SampleSemiMarkov.enable Adding key ", key, " interval ",
        	interval, " when ", firing_time, " dist ", clock)
		record = TransitionRecord(interval, handle)
		propagator.transition_state[key] = record
	end
    @debug("SampleSemiMarkov.enable exit")
end


# Remove a transition from the queue because it was disabled.
function Disable(propagator::NextReactionHazards, key, now,
        updated, rng)
	record = propagator.transition_state[key]
	# We store distributions in order to calculate remaining hazard
	# which will happen AFTER the state has changed.
	update!(propagator.firing_queue, record.heap_handle,
		OrderedSample(key, -1.))
	todelete = pop!(propagator.firing_queue)
	@assert(todelete.key == key && todelete.time == -1)
    if updated == :Disabled
    	record.heap_handle = -1 # This is the official sign it was disabled.
    elseif updated == :Fired
        # Deleting the key is slower for small, finite systems,
        # but it makes infinite (meaning long-running) systems possible.
        delete!(propagator.transition_state, key)
    else
        assert(updated == :Disabled || updated == :Fired)
    end
end


function print_next_reaction_hazards(propagator::NextReactionHazards)
    arr = Any[]
    for x in keys(propagator.transition_state)
        push!(arr, x)
    end
    sort!(arr)
    for trans in arr
        rec = propagator.transition_state[trans]
        if rec.distribution !== nothing
	        p=parameters(rec.distribution)
    	end
    end
end
