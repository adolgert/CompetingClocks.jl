# ----------------------------------------------------------------------
# Watchers
#
# The context carries zero or more *observer* middleware (likelihood
# trackers, debug recorders, ...) in a single heterogeneous `watchers`
# tuple. All fan-out is compile-time tuple recursion, so an absent
# watcher costs nothing and a present one is a fully-inlined call.
#
# `crn` (common random numbers) and `delayed` are NOT watchers: `crn`
# is a sampler decorator and `delayed` is per-clock state, so they stay
# dedicated fields.
#
# These helpers are internal; none are exported or public.
# ----------------------------------------------------------------------

# Watchers that want the FULL vector of distributions on a vector-enable
# (for importance sampling / likelihoods). Everything else receives the
# single sampled distribution, matching the sampler.
const LikelihoodWatcher = Union{TrackWatcher, TrajectoryWatcher, PathLikelihoods}

# --- fan-out: enable! -------------------------------------------------
# `dists` is what a LikelihoodWatcher sees (the full vector on a vector
# enable, the scalar otherwise); `sampled` is what everyone else sees.
@inline watch_enable_one!(w::LikelihoodWatcher, clock, dists, sampled, te, when, rng) =
    enable!(w, clock, dists, te, when, rng)
@inline watch_enable_one!(w, clock, dists, sampled, te, when, rng) =
    enable!(w, clock, sampled, te, when, rng)

@inline watch_enable!(::Tuple{}, clock, dists, sampled, te, when, rng) = nothing
@inline function watch_enable!(ws::Tuple, clock, dists, sampled, te, when, rng)
    watch_enable_one!(first(ws), clock, dists, sampled, te, when, rng)
    watch_enable!(Base.tail(ws), clock, dists, sampled, te, when, rng)
end

# --- fan-out: disable! / fire! / reset! / copy_clocks! ----------------
@inline watch_disable!(::Tuple{}, clock, when) = nothing
@inline function watch_disable!(ws::Tuple, clock, when)
    disable!(first(ws), clock, when)
    watch_disable!(Base.tail(ws), clock, when)
end

@inline watch_fire!(::Tuple{}, clock, when) = nothing
@inline function watch_fire!(ws::Tuple, clock, when)
    fire!(first(ws), clock, when)
    watch_fire!(Base.tail(ws), clock, when)
end

watch_reset!(::Tuple{}) = nothing
function watch_reset!(ws::Tuple)
    reset!(first(ws))
    watch_reset!(Base.tail(ws))
end

@inline watch_copy_clocks!(::Tuple{}, ::Tuple{}) = nothing
@inline function watch_copy_clocks!(dst::Tuple, src::Tuple)
    copy_clocks!(first(dst), first(src))
    watch_copy_clocks!(Base.tail(dst), Base.tail(src))
end

# --- compile-time lookup by watcher type ------------------------------
# Returns the first watcher that is an instance of `U`, or `nothing`.
# The `isa` checks are on concrete types, so they constant-fold and the
# result type is inferred exactly (the watcher type or `Nothing`).
@inline _first_of(::Type{U}, ::Tuple{}) where {U} = nothing
@inline function _first_of(::Type{U}, ws::Tuple) where {U}
    first(ws) isa U ? first(ws) : _first_of(U, Base.tail(ws))
end

# Filter a tuple, dropping `nothing` entries. Used at construction to
# assemble the watcher tuple from the optional middleware.
@inline _compact(::Tuple{}) = ()
@inline _compact(t::Tuple) =
    first(t) === nothing ? _compact(Base.tail(t)) : (first(t), _compact(Base.tail(t))...)


mutable struct SamplingContext{K,T,Sampler<:SSA,RNG,W<:Tuple,CRN,DS}
    sampler::Sampler # The actual sampler (may use internal key type)
    rng::RNG
    watchers::W      # Tuple of observer middleware (likelihood, debug, ...)
    crn::CRN
    split_weight::Float64
    time::T
    fixed_start::T
    sample_distribution::Int  # Given a vector of distributions during enabling, sample from this one.
    delayed::DS               # Either `Nothing` or `DelayedState{K,T}`
end


"""
    likelihood_of(ctx)

Return the likelihood watcher (`TrackWatcher`/`TrajectoryWatcher`/`PathLikelihoods`)
held by the context, or `nothing`. Compile-time tuple lookup; type-stable.
Internal accessor.
"""
@inline likelihood_of(ctx::SamplingContext) = _first_of(LikelihoodWatcher, ctx.watchers)

"""
    debug_of(ctx)

Return the `DebugWatcher` held by the context, or `nothing`. Internal accessor.
"""
@inline debug_of(ctx::SamplingContext) = _first_of(DebugWatcher, ctx.watchers)


"""
    internal_key(ctx, clock)

Map a user clock key to the sampler's internal key. Identity for a
non-delayed context (`DS === Nothing`); `(clock, :regular)` for a delayed
context. Dispatch-based, so it constant-folds.
"""
@inline internal_key(::SamplingContext{K,T,S,R,W,C,Nothing}, clock) where {K,T,S,R,W,C} = clock
@inline internal_key(::SamplingContext{K,T,S,R,W,C,DS}, clock) where {K,T,S,R,W,C,DS<:DelayedState} = (clock, :regular)


"""
    sampler_enable!(ctx, ikey, dist, te, when)

Enable `ikey` on the underlying sampler, wrapping the draw in common
random numbers when the context has a CRN recorder. This is the single
place the CRN branch lives; every enable path routes through here.
"""
@inline function sampler_enable!(ctx::SamplingContext, ikey, dist, te, when)
    crn = ctx.crn
    if crn !== nothing
        with_common_rng(crn, ikey, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, ikey, dist, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, ikey, dist, te, when, ctx.rng)
    end
end


"""
    SamplingContext(builder::SamplerBuilder, rng)

Uses the [`SamplerBuilder`](@ref) to make a SamplingContext.

`K` is always the *user* key type (`builder.clock_type`). In a delayed context
(`builder.support_delayed == true`) the public event identity is the tuple
`(clock, phase)`, with `phase` one of `:regular`, `:initiate`, or `:complete`.
The sampler and middleware store keys as `K_int`, which is `Tuple{K,Symbol}` for
delayed contexts and `K` for regular ones.
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

    # Observer middleware collapses into a single heterogeneous tuple.
    watchers = _compact((likelihood, debug))

    # Delayed state
    delayed_state = builder.support_delayed ? DelayedState{K,T}() : nothing

    SamplingContext{K,T,typeof(sampler),R,typeof(watchers),typeof(crn),typeof(delayed_state)}(
        sampler, rng, watchers, crn,
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
All subcomponents (sampler, watchers, CRN, delayed state) are cloned.
"""
function clone(sc::SamplingContext{K,T,Sampler,RNG,W,CRN,DS}, rng::RNG) where {K,T,Sampler,RNG,W,CRN,DS}
    delayed_copy = sc.delayed === nothing ? nothing : clone(sc.delayed)
    SamplingContext{K,T,Sampler,RNG,W,CRN,DS}(
        clone(sc.sampler),
        rng,
        map(clone, sc.watchers),
        isnothing(sc.crn) ? nothing : clone(sc.crn),
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
    likelihood = likelihood_of(ctx)
    if dist_index > 1 && !(likelihood isa PathLikelihoods)
        error("Can't sample from a later distribution unless likelihood_cnt>1 in SamplerBuilder")
    end
    dist_cnt = _likelihood_cnt(likelihood)
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
# enable! – regular clocks
#
# A single scalar body and a single vector body serve BOTH delayed and
# non-delayed contexts; `internal_key` folds the (clock, :regular) mapping.
# ----------------------------------------------------------------------

"""
    enable!(ctx, clock, dist, relative_te)

Enable a regular clock. In a delayed-support context the user key `clock`
is mapped to the internal key `(clock, :regular)`.
"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::UnivariateDistribution, relative_te::T) where {K,T}
    when = ctx.time
    te   = when + relative_te
    ikey = internal_key(ctx, clock)
    watch_enable!(ctx.watchers, ikey, dist, dist, te, when, ctx.rng)
    sampler_enable!(ctx, ikey, dist, te, when)
end

"""
    enable!(ctx, clock, dist)

Enable with no shift in enabling time.
"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::UnivariateDistribution) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end

"""
    enable!(ctx, clock, dist::Vector, relative_te)

Vectorized enable (importance sampling). The likelihood watcher receives the
full vector of distributions while the sampler and other watchers receive the
selected `dist[sample_idx]`.
"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector, relative_te::T) where {K,T}
    when = ctx.time
    te   = when + relative_te
    sample_idx = length(dist) == 1 ? 1 : ctx.sample_distribution
    sampled = dist[sample_idx]
    ikey = internal_key(ctx, clock)
    watch_enable!(ctx.watchers, ikey, dist, sampled, te, when, ctx.rng)
    sampler_enable!(ctx, ikey, sampled, te, when)
end

"""
    enable!(ctx, clock, dist::Vector)

Vectorized enable with no shift in enabling time.
"""
function enable!(ctx::SamplingContext{K,T}, clock::K, dist::Vector) where {K,T}
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
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS},
                 clock::K, delayed::Delayed, relative_te::T=zero(T)) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    # Clean up any existing events for this clock (e.g., if re-enabling during completion phase)
    if haskey(ctx.delayed.durations, clock)
        disable!(ctx, clock)
    end

    when = ctx.time
    te   = when + relative_te

    # Store duration distribution, not a sampled duration
    ctx.delayed.durations[clock] = delayed.duration

    ikey = (clock, :initiate)
    watch_enable!(ctx.watchers, ikey, delayed.initiation, delayed.initiation, te, when, ctx.rng)
    sampler_enable!(ctx, ikey, delayed.initiation, te, when)
end

"""
    enable!(ctx, clock, p::Pair, relative_te)

Enable a delayed reaction from a pair of distributions. Writing `d1 => d2`
builds an ordinary `Pair`, which this method converts to a `Delayed` and
forwards to the `Delayed` method above. This is how the user-facing syntax
`enable!(ctx, clock, Exponential(1.0) => Gamma(2.0))` is supported without any
type piracy on `Base.:(=>)`.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS},
                 clock::K, p::Pair{<:UnivariateDistribution,<:UnivariateDistribution},
                 relative_te::T=zero(T)) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    enable!(ctx, clock, Delayed(p), relative_te)
end

"""
    enable!(ctx, clock, p::Pair, relative_te)

A `Pair` of distributions specifies a delayed reaction, which this context was
not built to support. Guard method so the mistake fails at the API boundary
with a pointer to the fix, instead of deep inside a sampler.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,Nothing},
                 clock::K, p::Pair{<:UnivariateDistribution,<:UnivariateDistribution},
                 relative_te::T=zero(T)) where {K,T,Sampler,RNG,W,CRN}
    error("A Pair of distributions, initiation => duration, describes a delayed " *
          "reaction. Create the context with support_delayed=true to use it.")
end


# ----------------------------------------------------------------------
# disable!
# ----------------------------------------------------------------------

"""
    disable!(ctx, clock)

Non-delayed context: disable a single clock.
"""
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,Nothing}, clock::K) where {K,T,Sampler,RNG,W,CRN}
    when = ctx.time
    watch_disable!(ctx.watchers, clock, when)
    disable!(ctx.sampler, clock, when)
end

"""
    disable!(ctx, clock)

Delayed context: disable all phases associated with `clock` and clear any stored
duration distribution.
"""
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS}, clock::K) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    when = ctx.time
    for phase in (:regular, :initiate, :complete)
        internal_clock = (clock, phase)
        # Only disable if this phase is actually enabled
        if isenabled(ctx.sampler, internal_clock)
            watch_disable!(ctx.watchers, internal_clock, when)
            disable!(ctx.sampler, internal_clock, when)
        end
    end
    delete!(ctx.delayed.durations, clock)
end


# ----------------------------------------------------------------------
# next and next_delayed
# ----------------------------------------------------------------------

"""
    next(ctx::SamplingContext)

Return `(when, which)` for the next event from the underlying sampler, where
`which::K` is the user clock key.

The returned `(when, which)` is a *reservation*, not a commitment. It is valid
only until the next `enable!`, `disable!`, or `fire!` call changes the sampler's
state. Calling `next` twice without an intervening state change is undefined:
`CombinedNextReaction`'s `next` mutates internal sampler state. The two supported
responses are to fire it—call `fire!` with the returned clock and time—or to
decline it and stop the simulation (the fixed-horizon pattern). There is
deliberately no `peek`.

On a delayed context this method errors; use [`next_delayed`](@ref), which
returns `(when, clock, phase)`.
"""
function next(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,Nothing}) where {K,T,Sampler,RNG,W,CRN}
    next(ctx.sampler, ctx.time, ctx.rng)
end

function next(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS}) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    error("On a delayed context an event's identity is (clock, phase), so use " *
          "next_delayed(), which returns (when, clock, phase).")
end

"""
    next_delayed(ctx)

For delayed-support contexts, return `(when, which, phase)` where
`which::K` is the user key and `phase::Symbol` is `:regular`, `:initiate`,
or `:complete`.
"""
function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS}) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    when, internal_key = next(ctx.sampler, ctx.time, ctx.rng)
    clock, phase = internal_key
    return when, clock::K, phase::Symbol
end

function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,Nothing}) where {K,T,Sampler,RNG,W,CRN}
    error("next_delayed() requires support_delayed=true in the SamplerBuilder")
end


# ----------------------------------------------------------------------
# fire!
# ----------------------------------------------------------------------

"""
    fire!(ctx, clock, when)

Non-delayed context: fire a regular event and advance time.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,Nothing},
               clock::K, when::T) where {K,T,Sampler,RNG,W,CRN}
    watch_fire!(ctx.watchers, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.time = when
end

"""
    fire!(ctx, clock, when)

Delayed context: treat this as a regular (non-delayed) event with phase `:regular`.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS},
               clock::K, when::T) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    fire!(ctx, clock, :regular, when)
end

"""
    fire!(ctx, clock, phase, when)

Delayed context: fire an event with phase information.

- `phase == :regular`   : regular event
- `phase == :initiate`  : initiation of a delayed reaction; schedules completion
- `phase == :complete`  : completion of a delayed reaction; clears delayed state
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS},
               clock::K, phase::Symbol, when::T) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    internal_clock = (clock, phase)

    watch_fire!(ctx.watchers, internal_clock, when)
    fire!(ctx.sampler, internal_clock, when)

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
function enable_completion!(ctx::SamplingContext{K,T,Sampler,RNG,W,CRN,DS},
                            clock::K,
                            duration_dist::UnivariateDistribution,
                            when::T) where {K,T,Sampler,RNG,W,CRN,DS<:DelayedState{K,T}}
    ikey = (clock, :complete)
    te   = when   # duration distribution is measured from initiation time

    watch_enable!(ctx.watchers, ikey, duration_dist, duration_dist, te, when, ctx.rng)
    sampler_enable!(ctx, ikey, duration_dist, te, when)
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
    watch_reset!(ctx.watchers)
    reset!(ctx.sampler)
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
    watch_copy_clocks!(dst.watchers, src.watchers)
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

Return the set of enabled clock keys. In a delayed context each key is the
public event identity `(clock, phase)`, with `phase` one of `:regular`,
`:initiate`, or `:complete`—the same identity that `fire!` takes and
`next_delayed` returns.
"""
function enabled(ctx::SamplingContext)
    likelihood = likelihood_of(ctx)
    likelihood !== nothing && return enabled(likelihood)
    return enabled(ctx.sampler)
end


"""
    Base.length(sampling)

Total number of enabled clocks.
"""
function Base.length(ctx::SamplingContext)
    likelihood = likelihood_of(ctx)
    likelihood !== nothing && return length(likelihood)
    return length(ctx.sampler)
end


"""
    isenabled(sampling, clock)

Whether this user clock key is enabled. In delayed contexts this returns true
if any phase (`:regular`, `:initiate`, `:complete`) is enabled.
"""
function isenabled(ctx::SamplingContext{K}, clock::K) where {K}
    if ctx.delayed === nothing
        likelihood = likelihood_of(ctx)
        if likelihood !== nothing
            return isenabled(likelihood, clock)
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
In a delayed context, `which` is the public event identity `(clock, phase)`—the
same identity that `fire!` takes and `next_delayed` returns.
"""
function steploglikelihood(ctx::SamplingContext, when, which)
    likelihood = likelihood_of(ctx)
    if likelihood !== nothing
        return steploglikelihood(likelihood, ctx.time, when, which)
    elseif has_steploglikelihood(typeof(ctx.sampler))
        return steploglikelihood(ctx.sampler, ctx.time, when, which)
    else
        error("The sampler doesn't support steploglikelihood " *
              "unless you request it in the builder.")
    end
end


"""
    pathloglikelihood(ctx, endtime)

Path log-likelihood up to `endtime`, including the probability that no
event fires after the last event before `endtime`.
"""
function pathloglikelihood(ctx::SamplingContext, endtime)
    log_split = log(ctx.split_weight)
    likelihood = likelihood_of(ctx)
    if likelihood !== nothing
        @debug "Using likelihood object for trajectory"
        return pathloglikelihood(likelihood, endtime) .+ log_split
    elseif has_pathloglikelihood(typeof(ctx.sampler))
        return pathloglikelihood(ctx.sampler, endtime) .+ log_split
    else
        error("The sampler doesn't support pathloglikelihood " *
              "unless you request it in the builder.")
    end
end


"""
    enabled_history(ctx::SamplingContext)

Returns a `Vector{EnablingEntry{K,T}}` that has every time a clock was enabled.

See [`CompetingClocks.EnablingEntry`](@ref).
"""
function enabled_history(ctx::SamplingContext)
    debug = debug_of(ctx)
    if debug !== nothing
        return enabled_history(debug)
    else
        error("In order to get history of enabling create context with `recording=true`")
    end
end

"""
    disabled_history(ctx::SamplingContext)

Returns a `Vector{DisablingEntry{K,T}}` that has every time a clock was disabled.

See [`CompetingClocks.DisablingEntry`](@ref).
"""
function disabled_history(ctx::SamplingContext)
    debug = debug_of(ctx)
    if debug !== nothing
        return disabled_history(debug)
    else
        error("In order to get history of disabling create context with `recording=true`")
    end
end
