using DataStructures

struct MNRMTransition
    heap_handle::Int
    Lambda::Float64 # cumulative intensity
    distribution::UnivariateDistribution
    te::Float64  # Enabling time of distribution
    t0::Float64  # Enabling time of transition
end


"""
    MNRM by Anderson
"""
struct ModifiedNextReaction{T}
    # This stores everything we need to know about enabled transitions.
    # This struct requires the handle of the heap to always be greater than zero.
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    # This maps from transition to entry in the firing queue and remaining
    # integrated hazard. Has info on _disabled_ transitions.
    transition_entry::Dict{T,MNRMTransition}
end

function ModifiedNextReaction{T}() where {T}
    heap = MutableBinaryMinHeap{OrderedSample{T}}()
    ModifiedNextReaction{T}(heap, Dict{T,MNRMTransition}())
end