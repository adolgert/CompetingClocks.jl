export SamplingContext, enable!, fire!, isenabled, freeze_crn!
export sample_from_distribution!

"""
The SamplingContext is responsible for doing a perfect forwarding to multiple
components that implement the desired features in a sampler.

It uses a Context pattern where each type parameter is either present or
Nothing, and the compiler knows to optimize-away the Nothing.
We keep the internal logic simple. If something is present, call it.

 - memory
 - delayed transitions
 - hierarchical samplers
"""
mutable struct SamplingContext{K,T,Sampler<:SSA{K,T},RNG,Like,CRN,Dbg}
    sampler::Sampler # The actual sampler
    rng::RNG
    likelihood::Like
    crn::CRN
    debug::Dbg      # Union{Nothing, TrackingState{K,T}}
    split_weight::Float64
    time::T
    fixed_start::T
    sample_distribution::Int  # Given a vector of distributions during enabling, sample from this one.
end


function SamplingContext(builder::SamplerBuilder, rng::R) where {R<:AbstractRNG}
    K = builder.clock_type
    T = builder.time_type
    sampler = build_sampler(builder)
    if builder.likelihood_cnt > 1
        likelihood = PathLikelihoods{K,T}(builder.likelihood_cnt)
    elseif builder.trajectory_likelihood
        likelihood = TrajectoryWatcher{K,T}()
    else
        likelihood = nothing
    end
    if builder.step_likelihood && !has_steploglikelihood(typeof(sampler)) && isnothing(likelihood)
        likelihood = TrackWatcher{K,T}()
    end
    if builder.common_random
        crn = CommonRandom{K,R}()
    else
        crn = nothing
    end
    if builder.debug || builder.recording
        debug = DebugWatcher{K,T}(log=builder.debug)
    else
        debug = nothing
    end
    SamplingContext{K,T,typeof(sampler),R,typeof(likelihood),typeof(crn),typeof(debug)}(
        sampler, rng, likelihood, crn, debug, 1.0, builder.start_time, builder.start_time, 1
    )
end


function SamplingContext(::Type{K}, ::Type{T}, rng::R; kwargs...) where {K,T,R<:AbstractRNG}
    return SamplingContext(SamplerBuilder(K, T; kwargs...), rng)
end


"""
    clone(sampling, rng)

Given a `SamplingContext`, make a copy as though you had called the constructor
again. This is useful for creating a vector of `SamplingContext` for multi-threading
or for splitting paths. This clone will clone every part of the object.
Note that the random number generator is copied, and you will want that to be
unique for each clone.
    
You need to give each sampler its own random number generator, `rng`. If you
were making parallel RNGs with the `Random123` package, it might look like:

```julia
master_seed = (0xd897a239, 0x77ff9238)
rng = Philox4x((0, 0, 0, 0), master_seed)
Key = Int64
sampler = SamplingContext(SamplerBuilder(Key,Float64; trajectory_likelihood=true), rng)
observation_weight = zeros(Float64, particle_cnt)
total_weight = zeros(Float64, particle_cnt)
samplers = Vector{typeof(sampler)}(undef, particle_cnt)
for init_idx in 1:particle_cnt
    rng = Philox4x((0, 0, 0, init_idx), master_seed)
    samplers[init_idx] = clone(sampler, rng)
    init!(state[init_idx], samplers[init_idx]) # For some initialization of state.
end
````
"""
function clone(sc::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg}, rng::RNG) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg}(
        clone(sc.sampler),
        rng,
        isnothing(sc.likelihood) ? nothing : clone(sc.likelihood),
        isnothing(sc.crn) ? nothing : clone(sc.crn),
        isnothing(sc.debug) ? nothing : clone(sc.debug),
        1.0,
        sc.fixed_start,
        sc.fixed_start,
        sc.sample_distribution
    )
end


"""
    time(sampling)

Get the current simulation time. Simulation times are incremented by each
call to `fire!`.
"""
Base.time(ctx::SamplingContext) = ctx.time


"""
    sample_from_distribution!(sampling, index)

If you created a `SamplerBuilder` with `likelihood_cnt > 1`, then you are going
to call `enable!` for some or all clocks using a vector of distributions. In
that case, the `SamplingContext` will default to using the first distribution
to decide which event comes next. That way you get a sample from distribution 1
and likelihoods from distributions `(1, 2, and so on)`. For doing mixed
importance sampling, you will sample from several distributions and calculate
a weighted average. This functions lets you choose which of the enabled
distribution vector to sample.
"""
function sample_from_distribution!(ctx::SamplingContext, dist_index)
    if dist_index > 1 && !(ctx.likelihood isa PathLikelihoods)
        error("Can't sample from a later distribution unless likelihood_cnt>1 in SamplerBuilder")
    end
    dist_cnt = _likelihood_cnt(ctx.likelihood)
    if dist_index ∉ 1:dist_cnt
        error("Expected a distribution index between 1 and $dist_cnt")
    end
    ctx.sample_distribution = dist_index
end

"""
    freeze_crn!(ctx::SamplingContext)

After running the simulation to collect random draws, call this function to
stop collection of new draws and solely replay the draws that were collected,
using only fresh draws for clocks that weren't seen before.
"""
function freeze_crn!(ctx::SamplingContext)
    if ctx.crn !== nothing
        freeze_crn!(ctx.crn)
        ctx.time = ctx.fixed_start
    else
        error("Ask for common random numbers when creating the builder.")
    end
end


"""
    reset_crn!(ctx::SamplingContext)

This resets a sampler including erasing its stored common random numbers.
Using a simple `reset!(sampler)` won't erase the saved CRN draws.
"""
function reset_crn!(ctx::SamplingContext)
    if ctx.crn !== nothing
        reset_crn!(ctx.crn)
        ctx.time = ctx.fixed_start
    else
        error("Ask for common random numbers when creating the builder.")
    end
end


"""
    enable!(sampler, clock, distribution, relative_te)

Tell the sampler to start a clock with a time shift to the left or right.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Distributions.UnivariateDistribution`
 * `relative_te::TimeType` - The zero-point of the distribution relative to the current simulation
   time. A `-0.1` shifts the distribution to the left. A `0.1` says the clock
   cannot fire until after `0.1` units of time have passed.

"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist, relative_te::T) where {K,T}
    when = ctx.time
    te = when + relative_te
    ctx.likelihood !== nothing && enable!(ctx.likelihood, clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, clock, dist, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, clock, dist, te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, clock, dist, te, when, ctx.rng)
end


"""
    enable!(sampler, clock, distribution)

Tell the sampler to start a clock.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Distributions.UnivariateDistribution`


"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end


"""
    enable!(sampler, clock, distribution::Vector, relative_te)

Vectorized version of enable! for multiple-distributions in likelihood.
This will sample times from the first distribution in the vector unless
`sample_from_distribution!` is set. All distributions are used to calculate
a path likelihood for importance sampling.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Vector{Distributions.UnivariateDistribution}`
 * `relative_te::TimeType` - The zero-point of the distribution relative to the current simulation
   time. A `-0.1` shifts the distribution to the left. A `0.1` says the clock
   cannot fire until after `0.1` units of time have passed.

"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector, relative_te::T) where {K,T}
    when = ctx.time
    te = when + relative_te
    # This supports the ability of a caller to pass in a single distribution
    # for cases where they will not modify the importance of this event.
    sample_idx = length(dist) == 1 ? 1 : ctx.sample_distribution
    ctx.likelihood !== nothing && enable!(ctx.likelihood, clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, clock, dist[sample_idx], te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, clock, dist[sample_idx], te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, clock, dist[sample_idx], te, when, ctx.rng)
end



"""
    enable!(sampler, clock, distribution::Vector)

Vectorized version of enable! for multiple-distributions in likelihood.
This will sample times from the first distribution in the vector unless
`sample_from_distribution!` is set. All distributions are used to calculate
a path likelihood for importance sampling.

 * `sampler::SSA{KeyType,TimeType}` - The sampler to tell.
 * `clock::KeyType` - The ID of the clock. Can be a string, integer, tuple, etc.
 * `distribution::Vector{Distributions.UnivariateDistribution}`

This version enables the clock with no time-shifting.
"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end


"""
    disable!(sampler, clock)

Tell the sampler to forget a clock. The clock will be disabled at the time
set by the last call to `fire!`.
"""
function disable!(ctx::SamplingContext{K,T}, clock::K) where {K,T}
    when = ctx.time
    ctx.likelihood !== nothing && disable!(ctx.likelihood, clock, when)
    disable!(ctx.sampler, clock, when)
    ctx.debug !== nothing && disable!(ctx.debug, clock, when)
end


function next(ctx::SamplingContext{K,T}) where {K,T}
    next(ctx.sampler, ctx.time, ctx.rng)
end


"""
    fire!(ctx::SamplingContext, clock, when)

Decide the next `clock` key to fire and the time at which it fires. This sets
the internal time of the simulation to the new time `when`. It also disables
the given `clock` key.
"""
function fire!(ctx::SamplingContext{K,T}, clock::K, when::T) where {K,T}
    ctx.likelihood !== nothing && fire!(ctx.likelihood, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.debug !== nothing && fire!(ctx.debug, clock, when)
    ctx.time = when
end


"""
    reset!(sampling)

Clears all clock values in a `SamplingContext` at the top of a loop over
simulations.
"""
function reset!(ctx::SamplingContext{K,T}) where {K,T}
    ctx.likelihood !== nothing && reset!(ctx.likelihood)
    reset!(ctx.sampler)
    ctx.debug !== nothing && reset!(ctx.debug)
    ctx.time = ctx.fixed_start
end


"""
    copy_clocks!(dst::SamplingContext, src::SamplingContext)

Used for splitting a simulation near a rare event. You keep a vector of samplers
and once your current sampler reaches a simulation state, copy its state into
the vector of samplers and pick up where it left off.
"""
function copy_clocks!(dst::SamplingContext{K,T}, src::SamplingContext{K,T}) where {K,T}
    copy_clocks!(dst.sampler, src.sampler)
    jitter!(dst.sampler, ctx.when, ctx.rng)
    src.likelihood !== nothing && copy_clocks!(dst.likelihood, src.likelihood)
    src.debug !== nothing && copy_clocks!(dst.debug, src.debug)
    dst.split_weight = src.split_weight
    return dst
end


"""
    split!(dst::AbstractVector{S<:SamplingContext}, src::SamplingContext{K,T})

If the `src` sampler has gotten to a good spot in the simulation, split it into
multiple copies in the destination vector. Update the `split_weight` member of
each destination copy by `1/length(dst)`. The destination can be a view of
a vector that was created with `clone()`.
"""
function split!(dst::AbstractVector{S}, src::SamplingContext{K,T}) where {K,T,S<:SamplingContext}
    for cidx in eachindex(dst)
        copy_clocks(dst[cidx], src)
    end
    for sidx in eachindex(dst)
        dst[sidx].split_weight = src.split_weight / length(dst)
    end
end


"""
    enabled(sampling)

Returns a set of enabled clock keys. This will aggregate enabled clock keys
across hierarchical samplers.
"""
function enabled(ctx::SamplingContext)
    ctx.likelihood !== nothing && return enabled(ctx.likelihood)
    return enabled(ctx.sampler)
end


"""
    length(sampling)

The total number of enabled clocks.
"""
function Base.length(ctx::SamplingContext)
    ctx.likelihood !== nothing && return length(ctx.likelihood)
    return length(ctx.sampler)
end


"""
    isenabled(sampling, clock)

Boolean for whether this clock key is currently enabled, checked across
hierarchical samplers.
"""
function isenabled(ctx::SamplingContext{K}, clock::K) where {K}
    ctx.likelihood !== nothing && return isenabled(ctx.likelihood, clock)
    return isenabled(ctx.sampler, clock)
end

"""
    keytype(sampling)

The type for the clock key. If you can use a concrete type, that helps
performance.
"""
Base.keytype(ctx::SamplingContext{K}) where {K} = K

"""
    timetype(sampling)

The type for the time, usually a `Float64`.
"""
timetype(ctx::SamplingContext{K,T}) where {K,T} = T


"""
    steploglikelihood(sampling, when, which)

If the next event had clock key `which` and happened at time `when`,
what would be the log-likelihood of that event conditioned on the last event?
This calculates relative to the last call to `fire!()`. If you were to
integrate this value over all enabled events, from the last firing time to
`Inf` with the `QuadGK` package, it would sum to one.
"""
function steploglikelihood(ctx::SamplingContext, when, which)
    if ctx.likelihood !== nothing
        return steploglikelihood(ctx.likelihood, ctx.time, when, which)
    else
        try
            return steploglikelihood(ctx.sampler, ctx.time, when, which)
        catch
            error("The sampler doesn't support steploglikelihood " *
                  "unless you request it in the builder.")
        end
    end
end


"""
    pathloglikelihood(sampling, endtime)

Calculates the log-likelihood across all events since the start of the simulation
and up to the given end time. We include the end time because some statistics
are normalized not to the last event but to a specific time, so this includes
the probability that no event fired before the given end time.

This value can be used to weight paths for importance sampling.
"""
function pathloglikelihood(ctx::SamplingContext, endtime)
    log_split = log(ctx.split_weight)
    if ctx.likelihood !== nothing
        @debug "Using likelihood object for trajectory"
        return pathloglikelihood(ctx.likelihood, endtime) .+ log_split
    else
        try
            return pathloglikelihood(ctx.sampler, endtime) .+ log_split
        catch MethodError
            error("The sampler doesn't support pathloglikelihood " *
                  "unless you request it in the builder.")
        end
    end
end
