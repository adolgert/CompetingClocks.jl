# ----------------------------------------------------------------------
# Watchers
#
# The context carries zero or more *observer* middleware (likelihood
# trackers, debug recorders, ...) in a single heterogeneous `watchers`
# tuple. All fan-out is compile-time tuple recursion, so an absent
# watcher costs nothing and a present one is a fully-inlined call.
#
# `delayed` is NOT a watcher: it is per-clock state, so it stays a dedicated
# field. Common random numbers are no longer a context concern — the sampler owns
# its randomness as keyed streams (see keyed_streams.jl), and re-seeding /
# re-keying IS the common-random-numbers mechanism.
#
# These helpers are internal; none are exported or public.
# ----------------------------------------------------------------------

# Watchers that want the FULL vector of distributions on a vector-enable
# (for importance sampling / likelihoods). Everything else receives the
# single sampled distribution, matching the sampler.
const LikelihoodWatcher = Union{TrackWatcher, TrajectoryWatcher, PathLikelihoods}

# --- fan-out: enable! -------------------------------------------------
# `dists` is what a LikelihoodWatcher sees (the full vector on a vector
# enable, the scalar otherwise); `sampled` is what everyone else sees. No rng is
# threaded: watchers no longer take one (they never drew).
@inline watch_enable_one!(w::LikelihoodWatcher, clock, dists, sampled, te, when) =
    enable!(w, clock, dists, te, when)
@inline watch_enable_one!(w, clock, dists, sampled, te, when) =
    enable!(w, clock, sampled, te, when)

@inline watch_enable!(::Tuple{}, clock, dists, sampled, te, when) = nothing
@inline function watch_enable!(ws::Tuple, clock, dists, sampled, te, when)
    watch_enable_one!(first(ws), clock, dists, sampled, te, when)
    watch_enable!(Base.tail(ws), clock, dists, sampled, te, when)
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


mutable struct SamplingContext{K,T,Sampler<:SSA,RNG,W<:Tuple,DS}
    sampler::Sampler # The actual sampler (may use internal key type)
    # `rng` no longer threads into any sampler verb (the sampler owns keyed
    # streams). It is retained for context-owned, non-clock uses: it seeds the
    # sampler's streams at build time, and it draws a fresh seed for each split
    # copy in `split!`. It does not draw any clock's firing time.
    rng::RNG
    watchers::W      # Tuple of observer middleware (likelihood, debug, ...)
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
@inline internal_key(::SamplingContext{K,T,S,R,W,Nothing}, clock) where {K,T,S,R,W} = clock
@inline internal_key(::SamplingContext{K,T,S,R,W,DS}, clock) where {K,T,S,R,W,DS<:DelayedState} = (clock, :regular)


"""
    sampler_enable!(ctx, ikey, dist, te, when)

Enable `ikey` on the underlying sampler. The sampler owns its randomness (keyed
streams), so no rng is passed; every enable path routes through here.
"""
@inline function sampler_enable!(ctx::SamplingContext, ikey, dist, te, when)
    # Primal boundary: the sampler only draws numbers, so hand it a value-only
    # shadow of the distribution. Identity for ordinary Float64 parameters; an
    # AD extension strips tracer types so the sampler keeps running on Float64
    # while the differentiable copy lives on the likelihood watcher.
    dist = primal_distribution(dist)
    enable!(ctx.sampler, ikey, dist, te, when)
end


"""
    sampler_reenable!(ctx, ikey, dist, te, when)

Re-enable `ikey` on the underlying sampler, mirroring [`sampler_enable!`](@ref).
The sampler realizes the change under its own construction-time
[`coupling`](@ref); the context forwards no coupling. The sampler only draws
numbers, so it receives the primal (value-only) shadow of the distribution; the
differentiable copy stays on the likelihood watcher, which the caller fanned out
to separately.
"""
@inline function sampler_reenable!(ctx::SamplingContext, ikey, dist, te, when)
    dist = primal_distribution(dist)
    reenable!(ctx.sampler, ikey, dist, te, when)
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
    # The sampler owns its randomness. The context's rng is used ONCE here to
    # choose the sampler's stream seed, so a user still controls reproducibility
    # by the rng they pass, while the sampler thereafter draws from its own
    # per-clock streams rather than from this rng.
    rekey_streams!(sampler, rand(rng, UInt64))

    # Likelihood watcher. `L` is the accumulator number type: Float64 normally,
    # or a ForwardDiff.Dual eltype to carry derivatives through the likelihood.
    L = builder.likelihood_eltype
    if builder.likelihood_cnt > 1
        likelihood = PathLikelihoods{K_int,T,L}(builder.likelihood_cnt)
    elseif builder.path_likelihood
        likelihood = TrajectoryWatcher{K_int,T,L}()
    else
        likelihood = nothing
    end
    # Install a TrackWatcher for step likelihoods when either the sampler can't
    # compute steploglikelihood itself OR the caller asked for a non-Float64
    # accumulator. In the latter case the sampler's native steploglikelihood
    # would run on the PRIMAL distributions it received across the primal
    # boundary and silently return an underived Float64, dropping the very
    # derivatives the caller asked for. Forcing the watcher keeps the
    # Dual-parameterized distributions on the likelihood path.
    if builder.step_likelihood && isnothing(likelihood) &&
            (!has_steploglikelihood(typeof(sampler)) || L != Float64)
        likelihood = TrackWatcher{K_int,T}()
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

    SamplingContext{K,T,typeof(sampler),R,typeof(watchers),typeof(delayed_state)}(
        sampler, rng, watchers,
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

Clone a `SamplingContext` as though constructed again, using a new RNG. The
sampler and watchers are made fresh and empty (via [`similar_sampler`](@ref) and
the watchers' `clone`), and the new sampler's streams are re-keyed from `rng`, so
the clone runs an independent trajectory rather than a coupled copy. (For a
COUPLED copy of a live run, `clone` the sampler directly.)
"""
function clone(sc::SamplingContext{K,T,Sampler,RNG,W,DS}, rng::RNG) where {K,T,Sampler,RNG,W,DS}
    delayed_copy = sc.delayed === nothing ? nothing : clone(sc.delayed)
    fresh_sampler = similar_sampler(sc.sampler)
    rekey_streams!(fresh_sampler, rand(rng, UInt64))
    SamplingContext{K,T,Sampler,RNG,W,DS}(
        fresh_sampler,
        rng,
        map(clone, sc.watchers),
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


# Common random numbers are no longer a runtime mode toggled by freeze_crn! /
# reset_crn!. To couple two runs on the same randomness, build both contexts from
# the same seed (their samplers then draw identically per clock), or `clone` a
# sampler for a coupled continuation and `rekey_streams!` to decouple.


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
    watch_enable!(ctx.watchers, ikey, dist, dist, te, when)
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
    reenable!(ctx, clock, dist, relative_te)
    reenable!(ctx, clock, dist)

Re-evaluate the distribution of a clock that is CURRENTLY ENABLED, keeping its
age: the clock's remaining firing time becomes the law specified by the new
`(dist, te)` and the current time. This is the explicit form of calling
[`enable!`](@ref) again on an already-enabled clock. Which pathwise coupling
realizes the change — `:carry` (deterministic carry, the coupling pathwise/IPA
derivatives need) or `:redraw` (redraw-at-change) — is a property of the
underlying sampler, fixed at its construction and readable with
[`coupling`](@ref)`(ctx)`; the context forwards no coupling.

The watcher fan-out is IDENTICAL to `enable!`: the likelihood watcher and any
`TrajectoryRecorder` see a re-enable exactly as they see any enable, so they
close the OLD segment (accumulating its survival) and open the NEW one with the
new `(distribution, te)`. The piecewise-hazard likelihood and the recorder's
back-calculated firing uniform therefore use the CURRENT segment's distribution
automatically — the sampler's coupling affects only its retained draw, not the
recorded law.

`relative_te` re-anchors the enabling time exactly like `enable!`: the new
enabling time is `te = time(ctx) + relative_te`, so the clock's new age is
`-relative_te`. To KEEP the clock's age (the usual re-evaluation), pass the
clock's ORIGINAL enabling time as a shift, `relative_te = original_te −
time(ctx)` (a nonpositive number) — this is exactly what ChronoSim's
`sim_event_reenable` computes as `enable_time − sim.when`. The 3-argument form
mirrors `enable!(ctx, clock, dist)` and anchors at the current time
(`relative_te = 0`, age reset to 0), which is the shifted-re-anchoring case, not
the age-keeping one.
"""
function reenable!(ctx::SamplingContext{K,T}, clock::K, dist::UnivariateDistribution,
                   relative_te::T) where {K,T}
    when = ctx.time
    te   = when + relative_te
    ikey = internal_key(ctx, clock)
    watch_enable!(ctx.watchers, ikey, dist, dist, te, when)
    sampler_reenable!(ctx, ikey, dist, te, when)
end

function reenable!(ctx::SamplingContext{K,T}, clock::K, dist::UnivariateDistribution) where {K,T}
    reenable!(ctx, clock, dist, zero(T))
end

"""
    coupling(ctx::SamplingContext) -> Symbol

The re-evaluation coupling of the context's underlying sampler: `:carry` or
`:redraw`. This method forwards to the inner sampler (see [`coupling`](@ref))
so a recorder can read a run's coupling without reaching into the context.
"""
coupling(ctx::SamplingContext) = coupling(ctx.sampler)

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
    watch_enable!(ctx.watchers, ikey, dist, sampled, te, when)
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
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS},
                 clock::K, delayed::Delayed, relative_te::T=zero(T)) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    # Clean up any existing events for this clock (e.g., if re-enabling during completion phase)
    if haskey(ctx.delayed.durations, clock)
        disable!(ctx, clock)
    end

    when = ctx.time
    te   = when + relative_te

    # Store duration distribution, not a sampled duration
    ctx.delayed.durations[clock] = delayed.duration

    ikey = (clock, :initiate)
    watch_enable!(ctx.watchers, ikey, delayed.initiation, delayed.initiation, te, when)
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
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS},
                 clock::K, p::Pair{<:UnivariateDistribution,<:UnivariateDistribution},
                 relative_te::T=zero(T)) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    enable!(ctx, clock, Delayed(p), relative_te)
end

"""
    enable!(ctx, clock, p::Pair, relative_te)

A `Pair` of distributions specifies a delayed reaction, which this context was
not built to support. Guard method so the mistake fails at the API boundary
with a pointer to the fix, instead of deep inside a sampler.
"""
function enable!(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing},
                 clock::K, p::Pair{<:UnivariateDistribution,<:UnivariateDistribution},
                 relative_te::T=zero(T)) where {K,T,Sampler,RNG,W}
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
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}, clock::K) where {K,T,Sampler,RNG,W}
    when = ctx.time
    watch_disable!(ctx.watchers, clock, when)
    disable!(ctx.sampler, clock, when)
end

"""
    disable!(ctx, clock)

Delayed context: disable all phases associated with `clock` and clear any stored
duration distribution.
"""
function disable!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS}, clock::K) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
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
state. Calling `next` twice without an intervening state change returns valid
reservations, but they are not guaranteed to be identical across samplers:
`FirstReaction` redraws on each call, while `CombinedNextReaction` returns the
same cached reservation. Act on the most recent call's result. The two supported
responses are to fire it—call `fire!` with the returned clock and time—or to
decline it and stop the simulation (the fixed-horizon pattern). There is
deliberately no `peek`.

The underlying sampler requires a time argument that never decreases and never
advances past a pending event without firing it. The context satisfies this
automatically: it passes its own current time, which `fire!` keeps equal to the
last fired event's time.

On a delayed context this method errors; use [`next_delayed`](@ref), which
returns `(when, clock, phase)`.
"""
function next(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}) where {K,T,Sampler,RNG,W}
    next(ctx.sampler, ctx.time)
end

function next(ctx::SamplingContext{K,T,Sampler,RNG,W,DS}) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    error("On a delayed context an event's identity is (clock, phase), so use " *
          "next_delayed(), which returns (when, clock, phase).")
end

"""
    next_delayed(ctx)

For delayed-support contexts, return `(when, which, phase)` where
`which::K` is the user key and `phase::Symbol` is `:regular`, `:initiate`,
or `:complete`.
"""
function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,W,DS}) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    when, internal_key = next(ctx.sampler, ctx.time)
    clock, phase = internal_key
    return when, clock::K, phase::Symbol
end

function next_delayed(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}) where {K,T,Sampler,RNG,W}
    error("next_delayed() requires support_delayed=true in the SamplerBuilder")
end


# ----------------------------------------------------------------------
# fire!
# ----------------------------------------------------------------------

"""
    fire!(ctx, clock, when)

Non-delayed context: fire a regular event and advance time.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing},
               clock::K, when::T) where {K,T,Sampler,RNG,W}
    watch_fire!(ctx.watchers, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.time = when
end

"""
    fire!(ctx, clock, when)

Delayed context: treat this as a regular (non-delayed) event with phase `:regular`.
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS},
               clock::K, when::T) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    fire!(ctx, clock, :regular, when)
end

"""
    fire!(ctx, clock, phase, when)

Delayed context: fire an event with phase information.

- `phase == :regular`   : regular event
- `phase == :initiate`  : initiation of a delayed reaction; schedules completion
- `phase == :complete`  : completion of a delayed reaction; clears delayed state
"""
function fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS},
               clock::K, phase::Symbol, when::T) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
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
    force_fire!(ctx, clock, tstar)

Non-delayed context: fire the CHOSEN `clock` at the CHOSEN time `tstar`
regardless of the race, the branch step of the weak-derivative estimator, and
advance time to `tstar`. The watcher fan-out is the SAME as `fire!`, so an
attached `TrajectoryRecorder` records the forced firing with a well-defined
retained draw — for a forced firing the back-calculated uniform is still
`ccdf(dist, tstar - te)` — and cannot tell an imposed firing from a raced one.
The sampler's losers are re-conditioned on survival past `tstar` by
[`force_fire!`](@ref) on the underlying sampler.

The sampler's loser-redraw path draws from each loser's OWN keyed stream, so a
forced firing stays coupled across two contexts built from the same seed.
Requires a sampler with `supports_force == true`; otherwise the underlying verb
throws an `ArgumentError`.
"""
function force_fire!(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing},
                     clock::K, tstar::T) where {K,T,Sampler,RNG,W}
    watch_fire!(ctx.watchers, clock, tstar)
    force_fire!(ctx.sampler, clock, tstar)
    ctx.time = tstar
end

"""
    enabled_ages(ctx::SamplingContext, when=time(ctx)) -> Vector{Tuple{K,Float64}}

Every enabled clock paired with its age `when - te`, sorted by key — the
context-level form of the sampler verb, so an estimator taking its branching
decision at a context's decision point need not reach the raw sampler.
Defined for the non-delayed context, where the sampler's keys are the caller's
clock keys (`internal_key` is the identity); a delayed context wraps keys and
would need the reverse mapping before this can be offered there. Requires a
sampler with `supports_enabled_ages == true`.
"""
enabled_ages(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}) where {K,T,Sampler,RNG,W} =
    enabled_ages(ctx, time(ctx))
enabled_ages(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}, when::T) where {K,T,Sampler,RNG,W} =
    enabled_ages(ctx.sampler, when)

"""
    getindex(ctx::SamplingContext, clock) -> T

The scheduled (putative) firing time of the enabled `clock` — the context-level
form of the sampler's `getindex`, so an estimator that needs a clock's schedule
(the runner-up's residual in the smoothed-perturbation-analysis weight) need not
reach the raw sampler. Defined for the non-delayed context, where the sampler's
keys are the caller's clock keys. Throws `KeyError` when `clock` is not
currently enabled. Requires a scheduling sampler that stores putative times
(e.g. `CombinedNextReaction`); a rejection-style sampler has no schedule to
report.
"""
Base.getindex(ctx::SamplingContext{K,T,Sampler,RNG,W,Nothing}, clock::K) where {K,T,Sampler,RNG,W} =
    ctx.sampler[clock]

"""
    enable_completion!(ctx, clock, duration_dist, when)

Internal helper: schedule the completion event for a delayed reaction at time
`when + X` where `X` is drawn from `duration_dist`.
"""
function enable_completion!(ctx::SamplingContext{K,T,Sampler,RNG,W,DS},
                            clock::K,
                            duration_dist::UnivariateDistribution,
                            when::T) where {K,T,Sampler,RNG,W,DS<:DelayedState{K,T}}
    ikey = (clock, :complete)
    te   = when   # duration distribution is measured from initiation time

    watch_enable!(ctx.watchers, ikey, duration_dist, duration_dist, te, when)
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

Copy enabled clocks, stream state, and delayed state from `src` into `dst`. This
is a FAITHFUL copy: `dst`'s sampler inherits `src`'s keyed streams (generator
states and counts), so on its own it would reproduce `src`'s draws. To make the
copy diverge — the splitting use case — pair this with [`rekey_streams!`](@ref)
and `jitter!` as [`split!`](@ref) does.
"""
function copy_clocks!(dst::SamplingContext{K,T}, src::SamplingContext{K,T}) where {K,T}
    copy_clocks!(dst.sampler, src.sampler)
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
    jitter!(ctx::SamplingContext, when)

Resample every currently-enabled clock of the context's sampler from its own
keyed streams, forwarding to [`jitter!`](@ref) on the inner sampler. This pairs
with [`rekey_streams!`](@ref) to make a faithful copy diverge — exactly the
sequence [`split!`](@ref) performs internally — so branching code driving a
context need not reach into `ctx.sampler`. Jittering without a preceding re-key
reproduces the original draws.
"""
jitter!(ctx::SamplingContext{K,T}, when::T) where {K,T} = jitter!(ctx.sampler, when)


"""
    split!(dst, src)

Split a `src` context into multiple copies in `dst`, adjusting `split_weight`.
Each copy first inherits `src`'s clocks and streams via [`copy_clocks!`](@ref),
then is given its OWN randomness: `rekey_streams!` reseeds the copy's streams
(from that copy's context rng) and `jitter!` resamples the currently-enabled
clocks off the freshly-seeded streams. This replaces the old jitter-with-a-fresh-
rng differentiation; re-keying is the mechanism now that the sampler owns its
randomness.
"""
function split!(dst::AbstractVector{S}, src::SamplingContext{K,T}) where {K,T,S<:SamplingContext}
    for cidx in eachindex(dst)
        copy_clocks!(dst[cidx], src)
        rekey_streams!(dst[cidx].sampler, rand(dst[cidx].rng, UInt64))
        jitter!(dst[cidx].sampler, src.time)
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
