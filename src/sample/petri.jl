export Petri

"""
    Petri{KeyType,TimeType}()

If you want to test a simulation, it can be helpful to test unlikely events.
This sampler adopts the Petri net rule for which clock fires next: it's randomly
chosen among all enabled events. The returned time is always the previous time
plus one.
"""
mutable struct Petri{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    time_duration::Float64
    Petri{K,T}(dt=1.0) where {K,T} = new(Dict{K,EnablingEntry{K,T}}(), dt)
end

function Base.copy!(dst::Petri, src::Petri)
    copy!(dst.enabled, src.enabled)
    dst.time_duration = src.time_duration
    return dst
end

function next(propagator::Petri{K,T}, when::T, rng::AbstractRNG) where {K,T}
    chosen = rand(rng, keys(propagator.enabled))
    (when + propagator.time_duration, chosen)
end
