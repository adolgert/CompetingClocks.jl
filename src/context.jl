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


Base.time(ctx::SamplingContext) = ctx.time


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

function enable!(ctx::SamplingContext{K,T}, clock::K, dist) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end


# Vectorized version of enable! for multiple-distributions in likelihood
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector, relative_te::T) where {K,T}
    when = ctx.time
    te = when + relative_te
    ctx.likelihood !== nothing && enable!(ctx.likelihood, clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, clock, dist[ctx.sample_distribution], te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, clock, dist[ctx.sample_distribution], te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, clock, dist[ctx.sample_distribution], te, when, ctx.rng)
end


function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end


function disable!(ctx::SamplingContext{K,T}, clock::K) where {K,T}
    when = ctx.time
    ctx.likelihood !== nothing && disable!(ctx.likelihood, clock, when)
    disable!(ctx.sampler, clock, when)
    ctx.debug !== nothing && disable!(ctx.debug, clock, when)
end


function next(ctx::SamplingContext{K,T}) where {K,T}
    next(ctx.sampler, ctx.time, ctx.rng)
end


function fire!(ctx::SamplingContext{K,T}, clock::K, when::T) where {K,T}
    ctx.likelihood !== nothing && fire!(ctx.likelihood, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.debug !== nothing && fire!(ctx.debug, clock, when)
    ctx.time = when
end


function reset!(ctx::SamplingContext{K,T}) where {K,T}
    ctx.likelihood !== nothing && reset!(ctx.likelihood)
    reset!(ctx.sampler)
    ctx.debug !== nothing && reset!(ctx.debug)
    ctx.time = ctx.fixed_start
end


function Base.copy!(dst::SamplingContext{K,T}, src::SamplingContext{K,T}) where {K,T}
    copy!(dst.sampler, src.sampler)
    copy!(dst.rng, src.rng) # Will need to initialize this separately.
    src.likelihood !== nothing && copy!(dst.likelihood, src.likelihood)
    src.debug !== nothing && copy!(dst.debug, src.debug)
    dst.split_weight = src.split_weight
    return dst
end


function enabled(ctx::SamplingContext)
    ctx.likelihood !== nothing && return enabled(ctx.likelihood)
    return enabled(ctx.sampler)
end


function Base.length(ctx::SamplingContext)
    ctx.likelihood !== nothing && return length(ctx.likelihood)
    return length(ctx.sampler)
end


function isenabled(ctx::SamplingContext{K}, clock::K) where {K}
    ctx.likelihood !== nothing && return isenabled(ctx.likelihood, clock)
    return isenabled(ctx.sampler, clock)
end

Base.keytype(ctx::SamplingContext{K}) where {K} = K
timetype(ctx::SamplingContext{K,T}) where {K,T} = T


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


function trajectoryloglikelihood(ctx::SamplingContext, endtime)
    if ctx.likelihood !== nothing
        @debug "Using likelihood object for trajectory"
        return trajectoryloglikelihood(ctx.likelihood, endtime)
    else
        try
            return trajectoryloglikelihood(ctx.sampler, endtime)
        catch MethodError
            error("The sampler doesn't support trajectoryloglikelihood " *
                  "unless you request it in the builder.")
        end
    end
end
