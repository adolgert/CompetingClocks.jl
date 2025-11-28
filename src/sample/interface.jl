using Random: AbstractRNG
using Distributions: UnivariateDistribution
import Base: getindex, keys, length, keytype, haskey

export SSA, enable!, disable!, next, copy_clocks!,
    getindex, keys, length, keytype, reset!

"""
    SSA{KeyType,TimeType}

This abstract type represents a stochastic simulation algorithm (SSA). It is
parametrized by the clock ID, or key, and the type used for the time, which
is typically a Float64. The type of the key can be anything you would use
as a dictionary key. This excludes mutable values but includes a wide range
of identifiers useful for simulation. For instance, it could be a `String`,
but it could be a `Tuple{Int64,Int64,Int64}`, so that it indexes into a
complicated simulation state.
"""
abstract type SSA{Key,Time} end


"""
    enable!(sampler, clock, distribution, enablingtime, currenttime, RNG)

Tell the sampler to start a clock.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Distributions.UnivariateDistribution`
 * `enablingtime::TimeType` - The zero time for the clock's distribution, in absolute time. Usually equal to `when`.
 * `when::TimeType` - The current time of the simulation.
 * `rng::AbstractRNG` - A random number generator.

These times are **absolute** since the start of the simulation. The current time
should be `when`. If you want to shift the distribution so that this event cannot
happen for a little while then choose `enablingtime > when`. If you want to modify
the distribution by shifting it left, then choose `enablingtime < when`. Usually,
`enablingtime == when`. It is also possible to always use `enablingtime == when`
and use the `truncated()` function to modify distributions.
"""
function enable!(
    sampler::SSA{K,T},
    clock::K,
    distribution::UnivariateDistribution,
    te::T, # enabling time
    when::T, # current simulation time
    rng::AbstractRNG
) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end

"""
    reset!(sampler)

After a sampler is used for a simulation run, it has internal state. This
function resets that internal state to the initial value in preparation
for another sample run.
"""
function reset!(sampler::SSA{K,T}) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    copy_clocks!(destination_sampler, source_sampler)

This copies the state of the source sampler to the destination sampler, replacing
the current state of the destination sampler. This is useful for splitting
techniques where you make copies of a simulation and restart it with different
random number generators.
"""
function copy_clocks!(sampler::SSA{K,T}) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    clone(sampler)

Given an existing sampler, make a copy that has the same type and same
constructor options but has no data in it. Use this to initialize an array
of samplers.
"""
function clone(sampler::SSA{K,T}) where {K,T}
    error("Clone not implemented for $(typeof(sampler))")
end


"""
    jitter!(sampler, when, rng)

Takes a sampler that is at a statistically-valid stopping time (this is a technical
term) and resamples all clocks currently-enabled so that this copy of the sampler
will yield different results from another copy of the sampler despite its internal
cache of clock data.
"""
function jitter!(sampler::SSA{K,T}, when::T, rng::AbstractRNG) where {K,T}
    error("jitter! not implemented for $(typeof(sampler))")
end


"""
    disable!(sampler, clock, when)

Tell the sampler to forget a clock. We include the current simulation time
because some Next Reaction methods use this to optimize sampling.
"""
function disable!(sampler::SSA{K,T}, clock::K, when::T) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end

"""
    next(sampler, when, rng)

Ask the sampler for what happens next, in the form of
`(when, which)::Tuple{TimeType,KeyType}`. `rng` is a random number generator.
"""
function next(sampler::SSA{K,T}, when::T, rng::AbstractRNG) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    getindex(sampler, clock::KeyType)

Return stored state for a particular clock. If the clock does not exist,
a `KeyError` will be thrown. Different samplers have **different stored state.**
"""
function Base.getindex(sampler::SSA{K,T}, clock::K) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    keys(sampler)

Return all stored clocks as a vector.
"""
function Base.keys(sampler::SSA)
    error("Not implemented for $(typeof(sampler))")
end


"""
    length(sampler)::Int64

Return the number of stored clocks.
"""
function Base.length(sampler::SSA)
    error("Not implemented for $(typeof(sampler))")
end


"""
    keytype(sampler)

Return the type of clock keys.
"""
Base.keytype(::SSA{K,T}) where {K,T} = K


"""
    timetype(sampler)

Return the type of clock times.
"""
timetype(::SSA{K,T}) where {K,T} = T


"""
    haskey(sampler, key)

Return a boolean.
"""
function Base.haskey(sampler::SSA{K,T}, key) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end


"""
    enabled(sampler)

Returns a read-only set of currently-enabled clocks.
"""
function enabled(sampler::SSA{K,T}) where {K,T}
    error("Not implemented for $(typeof(sampler))")
end
