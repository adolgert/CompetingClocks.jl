using DataStructures: MutableBinaryMinHeap, extract_all!, update!

export FirstToFire

# Helper struct for FirstToFire to be able to resample.
struct FTFEntry{T,D}
    handle::Int
    distribution::D
    te::T
end


"""
    FirstToFire{KeyType,TimeType}()

This sampler is often the fastest for non-exponential distributions.
When a clock is first enabled, this sampler asks the clock when it would
fire and saves that time in a sorted heap of future times. Then it works
through the heap, one by one. When a clock is disabled, its future firing time
is removed from the list. There is no memory of previous firing times.
"""
mutable struct FirstToFire{K,T} <: SSA{K,T}
    firing_queue::MutableBinaryMinHeap{OrderedSample{K,T}}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{K,FTFEntry{T,UnivariateDistribution}}
end


function FirstToFire{K,T}() where {K,T}
    heap = MutableBinaryMinHeap{OrderedSample{K,T}}()
    state = Dict{K,FTFEntry{T,UnivariateDistribution}}()
    FirstToFire{K,T}(heap, state)
end


clone(::FirstToFire{K,T}) where {K,T} = FirstToFire{K,T}()


function reset!(propagator::FirstToFire{K,T}) where {K,T}
    empty!(propagator.firing_queue)
    empty!(propagator.transition_entry)
end

function copy_clocks!(dst::FirstToFire{K,T}, src::FirstToFire{K,T}) where {K,T}
    dst.firing_queue = deepcopy(src.firing_queue)
    copy!(dst.transition_entry, src.transition_entry)
    return dst
end


function jitter!(propagator::FirstToFire{K,T}, when::T, rng::AbstractRNG) where {K,T}
    for (clock, entry) in propagator.transition_entry
        te = entry.te
        distribution = entry.distribution
        if te < when
            when_fire = te + rand(rng, truncated(distribution, when - te, typemax(T)))
        else
            when_fire = te + rand(rng, distribution)
        end
        update!(propagator.firing_queue, entry.handle, OrderedSample{K,T}(clock, when_fire))
    end
end


# Finds the next one without removing it from the queue.
function next(propagator::FirstToFire{K,T}, when::T, rng::AbstractRNG) where {K,T}
    least = if !isempty(propagator.firing_queue)
        first(propagator.firing_queue)
    else
        OrderedSample(nothing, typemax(T))
    end
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
        heap_handle = propagator.transition_entry[clock].handle
        update!(propagator.firing_queue, heap_handle, OrderedSample{K,T}(clock, when_fire))
    else
        heap_handle = push!(propagator.firing_queue, OrderedSample{K,T}(clock, when_fire))
    end
    propagator.transition_entry[clock] = FTFEntry{T,UnivariateDistribution}(
        heap_handle,
        distribution,
        te
    )
end


fire!(propagator::FirstToFire{K,T}, clock::K, when::T) where {K,T} = disable!(propagator, clock, when)


function disable!(propagator::FirstToFire{K,T}, clock::K, when::T) where {K,T}
    entry = get(propagator.transition_entry, clock, nothing)
    entry === nothing && throw(KeyError(clock))
    delete!(propagator.firing_queue, entry.handle)
    delete!(propagator.transition_entry, clock)
end

"""
    getindex(sampler::FirstToFire{K,T}, clock::K)

For the `FirstToFire` sampler, returns the stored firing time associated to the clock.
"""
function Base.getindex(propagator::FirstToFire{K,T}, clock::K) where {K,T}
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock].handle
        return getfield(propagator.firing_queue[heap_handle], :time)
    else
        throw(KeyError(clock))
    end
end

function Base.keys(propagator::FirstToFire)
    return collect(keys(propagator.transition_entry))
end

function Base.length(propagator::FirstToFire)
    return length(propagator.transition_entry)
end


Base.haskey(propagator::FirstToFire{K,T}, clock::K) where {K,T} = haskey(propagator.transition_entry, clock)
Base.haskey(propagator::FirstToFire{K,T}, clock) where {K,T} = false

enabled(propagator::FirstToFire) = keys(propagator.transition_entry)
isenabled(propagator::FirstToFire{K,T}, clock::K) where {K,T} = haskey(propagator.transition_entry, clock)
