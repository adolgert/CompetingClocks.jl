export Petri

"""
    Petri{KeyType,TimeType}()

If you want to test a simulation, it can be helpful to test unlikely events.
This sampler adopts the Petri net rule for which clock fires next: it's randomly
chosen among all enabled events. The returned time is always the previous time
plus one.
"""
mutable struct Petri{K,T} <: SSA{K,T}
    watcher::TrackWatcher{K,T}
    time_duration::Float64
    Petri{K,T}() where {K,T} = new(TrackWatcher{K,T}(), 1.0)
end

reset!(propagator::Petri{K,T}) where {K,T} = reset!(propagator.watcher)

Base.copy!(dst::Petri{K,T}, src::Petri{K,T}) where {K,T} = copy!(dst.watcher, src.watcher)

# Finds the next one without removing it from the queue.
function next(propagator::Petri{K,T}, when::T, rng::AbstractRNG) where {K,T}
    chosen = rand(rng, keys(propagator.watcher.enabled))
    (when + propagator.time_duration, chosen)
end


function enable!(
    propagator::Petri{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T, rng::AbstractRNG) where {K,T}
    return enable!(propagator.watcher, clock, distribution, te, when, rng)
end


function disable!(propagator::Petri{K,T}, clock::K, when::T) where {K,T}
    return disable!(propagator.watcher, clock, when)
end

"""
    getindex(sampler::Petri{K,T}, clock::K)

For the `Petri` sampler, returns the `EnablingEntry` time associated to the clock.
"""
function Base.getindex(propagator::Petri{K,T}, clock::K) where {K,T}
    return propagator.watcher.enabled[clock]
end

function Base.keys(propagator::Petri)
    return collect(keys(propagator.watcher.enabled))
end

function Base.length(propagator::Petri)
    return length(propagator.watcher)
end


Base.haskey(propagator::Petri, clock) = isenabled(propagator.watcher, clock)
