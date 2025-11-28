export SamplingContext, enable!, fire!, isenabled, freeze_crn!
export sample_from_distribution!
export next_delayed, timetype

mutable struct SamplingContext{K,T,Sampler<:SSA,RNG,Like,CRN,Dbg,DS}
    sampler::Sampler # The actual sampler (may use internal key type)
    rng::RNG
    likelihood::Like
    crn::CRN
    debug::Dbg      # Union{Nothing, TrackingState{K_int,T}}
    split_weight::Float64
    time::T
    fixed_start::T
    sample_distribution::Int  # Given a vector of distributions during enabling, sample from this one.
    delayed::DS               # Either `Nothing` or `DelayedState{K,T}`
end


"""
    SamplingContext(builder::SamplerBuilder, rng)

Uses the [`SamplerBuilder`](@ref) to make a SamplingContext.

`K` is always the *user* key type (`builder.clock_type`). The sampler and
middleware use an internal key type `K_int` which is equal to `K` for regular
contexts and `Tuple{K,Symbol}` when `builder.support_delayed == true`.
"""
function SamplingContext(builder::SamplerBuilder, rng::R) where {R<:AbstractRNG}
    K = builder.clock_type
    T = builder.time_type

    K_int = builder.support_delayed ? Tuple{K,Symbol} : K

    sampler = build_sampler(builder)

    # Likelihood watcher
    if builder.likelihood_cnt > 1
        likelihood = PathLikelihoods{K_int,T}(builder.likelihood_cnt)
    elseif builder.path_likelihood
        likelihood = TrajectoryWatcher{K_int,T}()
    else
        likelihood = nothing
    end
    if builder.step_likelihood && !has_steploglikelihood(typeof(sampler)) && isnothing(likelihood)
        likelihood = TrackWatcher{K_int,T}()
    end

    # Common random numbers
    if builder.common_random
        crn = CommonRandom{K_int,R}()
    else
        crn = nothing
    end

    # Debug / recording
    if builder.debug || builder.recording
        debug = DebugWatcher{K_int,T}(log=builder.debug)
    else
        debug = nothing
    end

    # Delayed state
    delayed_state = builder.support_delayed ? DelayedState{K,T}() : nothing

    SamplingContext{K,T,typeof(sampler),R,typeof(likelihood),typeof(crn),typeof(debug),typeof(delayed_state)}(
        sampler, rng, likelihood, crn, debug,
        1.0, builder.start_time, builder.start_time, 1,
        delayed_state,
    )
end


"""
    SamplingContext(::Type{K}, ::Type{T}, rng::AbstractRNG; kwargs...)

Convenience wrapper that allocates a `SamplerBuilder` then builds
a `SamplingContext`.
"""
function SamplingContext(::Type{K}, ::Type{T}, rng::R; kwargs...) where {K,T,R<:AbstractRNG}
    return SamplingContext(SamplerBuilder(K, T; kwargs...), rng)
end


"""
    clone(sampling, rng)

Clone a `SamplingContext` as though constructed again, using a new RNG.
All subcomponents (sampler, likelihood, CRN, debug, delayed state) are cloned.
"""
function clone(sc::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS}, rng::RNG) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS}
    delayed_copy = sc.delayed === nothing ? nothing : clone(sc.delayed)
    SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS}(
        clone(sc.sampler),
        rng,
        isnothing(sc.likelihood) ? nothing : clone(sc.likelihood),
        isnothing(sc.crn) ? nothing : clone(sc.crn),
        isnothing(sc.debug) ? nothing : clone(sc.debug),
        1.0,
        sc.fixed_start,
        sc.fixed_start,
        sc.sample_distribution,
        delayed_copy,
    )
end


"""
    time(sampling)

Current simulation time.
"""
Base.time(ctx::SamplingContext) = ctx.time


"""
    sample_from_distribution!(sampling, index)

Choose which element of a vector-of-distributions to sample from when using
path likelihoods (importance sampling).
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

Switch to using only previously recorded common random numbers.
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

Reset stored CRN draws and time.
"""
function reset_crn!(ctx::SamplingContext)
    if ctx.crn !== nothing
        reset_crn!(ctx.crn)
        ctx.time = ctx.fixed_start
    else
        error("Ask for common random numbers when creating the builder.")
    end
end


# ----------------------------------------------------------------------
# enable! – regular clocks (non-delayed contexts)
# ----------------------------------------------------------------------

"""
    enable!(ctx, clock, dist, relative_te)

Enable a regular clock in a non-delayed context.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing},
                 clock::K, dist, relative_te::T) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    when = ctx.time
    te   = when + relative_te
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

function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing},
                 clock::K, dist) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    enable!(ctx, clock, dist, zero(T))
end


"""
    enable!(ctx, clock, dist::Vector, relative_te)

Vectorized enable in a non-delayed context.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing},
                 clock::K, dist::Vector, relative_te::T) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    when = ctx.time
    te   = when + relative_te
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

function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing},
                 clock::K, dist::Vector) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    enable!(ctx, clock, dist, zero(T))
end


# ----------------------------------------------------------------------
# enable! – regular clocks (delayed contexts)
# ----------------------------------------------------------------------

"""
    enable!(ctx, clock, dist, relative_te)

Regular (non-delayed) clock enable in a delayed-support context.
User key `clock::K` is mapped to internal `(clock, :regular)`.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                 clock::K, dist, relative_te::T) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    when = ctx.time
    te   = when + relative_te
    internal_clock = (clock, :regular)

    ctx.likelihood !== nothing && enable!(ctx.likelihood, internal_clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, internal_clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, internal_clock, dist, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, internal_clock, dist, te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, internal_clock, dist, te, when, ctx.rng)
end

function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                 clock::K, dist) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    enable!(ctx, clock, dist, zero(T))
end


"""
    enable!(ctx, clock, dist::Vector, relative_te)

Vectorized enable for regular clocks in a delayed context.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                 clock::K, dist::Vector, relative_te::T) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    when = ctx.time
    te   = when + relative_te
    sample_idx = length(dist) == 1 ? 1 : ctx.sample_distribution
    internal_clock = (clock, :regular)

    ctx.likelihood !== nothing && enable!(ctx.likelihood, internal_clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, internal_clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, internal_clock, dist[sample_idx], te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, internal_clock, dist[sample_idx], te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, internal_clock, dist[sample_idx], te, when, ctx.rng)
end

function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                 clock::K, dist::Vector) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    enable!(ctx, clock, dist, zero(T))
end


# ----------------------------------------------------------------------
# enable! – delayed clocks
# ----------------------------------------------------------------------

"""
    enable!(ctx, clock, delayed::Delayed, relative_te)

Enable a delayed reaction. This:

  * Stores the **duration distribution** `delayed.duration` in
    `ctx.delayed.durations[clock]`.
  * Enables an initiation event with internal key `(clock, :initiate)`
    and distribution `delayed.initiation`.

The `relative_te` shift applies only to the initiation, not the completion.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                 clock::K, delayed::Delayed, relative_te::T=zero(T)) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    # Clean up any existing events for this clock (e.g., if re-enabling during completion phase)
    if haskey(ctx.delayed.durations, clock)
        disable!(ctx, clock)
    end

    when = ctx.time
    te   = when + relative_te

    # Store duration distribution, not a sampled duration
    ctx.delayed.durations[clock] = delayed.duration

    internal_clock = (clock, :initiate)

    ctx.likelihood !== nothing && enable!(ctx.likelihood, internal_clock, delayed.initiation, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, internal_clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, internal_clock, delayed.initiation, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, internal_clock, delayed.initiation, te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, internal_clock, delayed.initiation, te, when, ctx.rng)
end


# ----------------------------------------------------------------------
# disable!
# ----------------------------------------------------------------------

"""
    disable!(ctx, clock)

Non-delayed context: disable a single clock.
"""
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing}, clock::K) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    when = ctx.time
    ctx.likelihood !== nothing && disable!(ctx.likelihood, clock, when)
    disable!(ctx.sampler, clock, when)
    ctx.debug !== nothing && disable!(ctx.debug, clock, when)
end

"""
    disable!(ctx, clock)

Delayed context: disable all phases associated with `clock` and clear any stored
duration distribution.
"""
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS}, clock::K) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    when = ctx.time
    for phase in (:regular, :initiate, :complete)
        internal_clock = (clock, phase)
        # Only disable if this phase is actually enabled
        if isenabled(ctx.sampler, internal_clock)
            ctx.likelihood !== nothing && disable!(ctx.likelihood, internal_clock, when)
            disable!(ctx.sampler, internal_clock, when)
            ctx.debug !== nothing && disable!(ctx.debug, internal_clock, when)
        end
    end
    if ctx.delayed !== nothing
        delete!(ctx.delayed.durations, clock)
    end
end


# ----------------------------------------------------------------------
# next and next_delayed
# ----------------------------------------------------------------------

"""
    next(ctx::SamplingContext)

Return `(when, internal_key)` for the next event from the underlying sampler.
For non-delayed contexts, `internal_key` is the user key. For delayed contexts
it is the internal key, e.g. `(user_key, phase)`.
"""
function next(ctx::SamplingContext{K,T}) where {K,T}
    next(ctx.sampler, ctx.time, ctx.rng)
end

"""
    next_delayed(ctx)

For delayed-support contexts, return `(when, which, phase)` where
`which::K` is the user key and `phase::Symbol` is `:regular`, `:initiate`,
or `:complete`.
"""
function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS}) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    when, internal_key = next(ctx.sampler, ctx.time, ctx.rng)
    clock, phase = internal_key
    return when, clock::K, phase::Symbol
end

function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing}) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    error("next_delayed() requires support_delayed=true in the SamplerBuilder")
end


# ----------------------------------------------------------------------
# fire!
# ----------------------------------------------------------------------

"""
    fire!(ctx, clock, when)

Non-delayed context: fire a regular event and advance time.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,Nothing},
               clock::K, when::T) where {K,T,Sampler,RNG,Like,CRN,Dbg}
    ctx.likelihood !== nothing && fire!(ctx.likelihood, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.debug !== nothing && fire!(ctx.debug, clock, when)
    ctx.time = when
end

"""
    fire!(ctx, clock, when)

Delayed context: treat this as a regular (non-delayed) event with phase `:regular`.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
               clock::K, when::T) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    fire!(ctx, clock, :regular, when)
end

"""
    fire!(ctx, clock, phase, when)

Delayed context: fire an event with phase information.

- `phase == :regular`   : regular event
- `phase == :initiate`  : initiation of a delayed reaction; schedules completion
- `phase == :complete`  : completion of a delayed reaction; clears delayed state
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
               clock::K, phase::Symbol, when::T) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    internal_clock = (clock, phase)

    ctx.likelihood !== nothing && fire!(ctx.likelihood, internal_clock, when)
    fire!(ctx.sampler, internal_clock, when)
    ctx.debug !== nothing && fire!(ctx.debug, internal_clock, when)

    if phase === :initiate
        # Look up duration distribution and enable completion
        duration_dist = ctx.delayed.durations[clock]
        enable_completion!(ctx, clock, duration_dist, when)
    elseif phase === :complete
        # Completion finished; clear stored distribution
        delete!(ctx.delayed.durations, clock)
    end

    ctx.time = when
end

"""
    enable_completion!(ctx, clock, duration_dist, when)

Internal helper: schedule the completion event for a delayed reaction at time
`when + X` where `X` is drawn from `duration_dist`.
"""
function enable_completion!(ctx::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg,DS},
                            clock::K,
                            duration_dist::UnivariateDistribution,
                            when::T) where {K,T,Sampler,RNG,Like,CRN,Dbg,DS<:DelayedState{K,T}}
    internal_clock = (clock, :complete)
    dist = duration_dist
    te   = when   # duration distribution is measured from initiation time

    ctx.likelihood !== nothing && enable!(ctx.likelihood, internal_clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, internal_clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, internal_clock, dist, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, internal_clock, dist, te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, internal_clock, dist, te, when, ctx.rng)
end


# ----------------------------------------------------------------------
# reset!, copy_clocks!, split!, enabled, length, isenabled
# ----------------------------------------------------------------------

"""
    reset!(sampling)

Clear all enabled clocks (and delayed state if present) and reset time to
the fixed start time.
"""
function reset!(ctx::SamplingContext{K,T}) where {K,T}
    ctx.likelihood !== nothing && reset!(ctx.likelihood)
    reset!(ctx.sampler)
    ctx.debug !== nothing && reset!(ctx.debug)
    if ctx.delayed !== nothing
        reset!(ctx.delayed)
    end
    ctx.time = ctx.fixed_start
end


"""
    copy_clocks!(dst, src)

Copy enabled clocks and delayed state from `src` into `dst`, jittering the
destination sampler so it will generate different `next()` samples.
"""
function copy_clocks!(dst::SamplingContext{K,T}, src::SamplingContext{K,T}) where {K,T}
    copy_clocks!(dst.sampler, src.sampler)
    jitter!(dst.sampler, src.time, dst.rng)
    src.likelihood !== nothing && copy_clocks!(dst.likelihood, src.likelihood)
    src.debug      !== nothing && copy_clocks!(dst.debug, src.debug)
    dst.split_weight = src.split_weight

    if (src.delayed !== nothing) && (dst.delayed !== nothing)
        empty!(dst.delayed.durations)
        for (k,v) in src.delayed.durations
            dst.delayed.durations[k] = v
        end
    end

    return dst
end


"""
    split!(dst, src)

Split a `src` context into multiple copies in `dst`, adjusting `split_weight`.
"""
function split!(dst::AbstractVector{S}, src::SamplingContext{K,T}) where {K,T,S<:SamplingContext}
    for cidx in eachindex(dst)
        copy_clocks!(dst[cidx], src)
    end
    for sidx in eachindex(dst)
        dst[sidx].split_weight = src.split_weight / length(dst)
    end
end


"""
    enabled(sampling)

Return the set of enabled clock keys. For delayed contexts this returns the
*internal* keys, e.g. `(K, Symbol)`.
"""
function enabled(ctx::SamplingContext)
    ctx.likelihood !== nothing && return enabled(ctx.likelihood)
    return enabled(ctx.sampler)
end


"""
    Base.length(sampling)

Total number of enabled clocks.
"""
function Base.length(ctx::SamplingContext)
    ctx.likelihood !== nothing && return length(ctx.likelihood)
    return length(ctx.sampler)
end


"""
    isenabled(sampling, clock)

Whether this user clock key is enabled. In delayed contexts this returns true
if any phase (`:regular`, `:initiate`, `:complete`) is enabled.
"""
function isenabled(ctx::SamplingContext{K}, clock::K) where {K}
    if ctx.delayed === nothing
        if ctx.likelihood !== nothing
            return isenabled(ctx.likelihood, clock)
        else
            return isenabled(ctx.sampler, clock)
        end
    else
        return isenabled(ctx.sampler, (clock, :regular)) ||
               isenabled(ctx.sampler, (clock, :initiate)) ||
               isenabled(ctx.sampler, (clock, :complete))
    end
end


"""
    keytype(sampling)

User clock key type `K`.
"""
Base.keytype(ctx::SamplingContext{K}) where {K} = K

"""
    timetype(sampling)

Time type `T`, usually `Float64`.
"""
timetype(ctx::SamplingContext{K,T}) where {K,T} = T


"""
    steploglikelihood(ctx, when, which)

Step log-likelihood of an event `which` at `when`.
For delayed contexts, `which` should be the internal key (e.g. `(clock, phase)`).
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
    pathloglikelihood(ctx, endtime)

Path log-likelihood up to `endtime`, including the probability that no
event fires after the last event before `endtime`.
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
