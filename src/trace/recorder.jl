using Distributions: UnivariateDistribution, logccdf

# ---------------------------------------------------------------------------
# The recording watcher.
#
# A `TrajectoryRecorder` is a watcher (an `EnabledWatcher`, so it receives the
# same enable!/disable!/fire! verbs the samplers do) whose product is the
# ORDERED SEQUENCE OF FIRINGS of a trajectory. Each firing is stored together
# with the survival-space uniform `u` that satisfies the retained-draw
# identity
#
#     when == te + invlogccdf(dist, log(u))    equivalently    u == ccdf(dist, when - te).
#
# The recorder is deliberately SAMPLER-AGNOSTIC: it never asks the sampler for
# `u`. It knows `(dist, te)` from the enable! it received, and at fire time it
# BACK-CALCULATES `u` from the observed firing time. This is the same trick the
# next-reaction samplers use internally to recover "what the uniform would have
# been", promoted here to a recorded contract obligation. A first-reaction
# sampler that never draws by inversion and a next-reaction sampler that stores
# a survival uniform therefore produce byte-identical records for the same
# trajectory — the recorded firing sequence is a sufficient statistic for the
# path likelihood, so the score estimator and offline IPA read it and never
# touch sampler internals.
# ---------------------------------------------------------------------------

"""
    ClockFiredRecord{K,T}(clock, when, te, distribution, u, logu)

One firing as the recorder captured it: clock `clock` fired at absolute time
`when`, its lifetime distribution was `distribution` measured from enabling
time `te`, and `u` is the SURVIVAL-SPACE uniform of the clock's TOTAL lifetime
satisfying the retained-draw identity

```
when == te + invlogccdf(distribution, log(u))    equivalently    u == ccdf(distribution, when - te).
```

`logu = log(u)` is stored alongside `u` because the identity lives in
log-survival space: a deep tail underflows `u` to `0.0` while `logu` stays
finite, and it is `logu` (not `u`) that inverts back to `when` without loss.

`u` (and `logu`) are `NaN` when the enabling information needed to
back-calculate them was unavailable — this should not happen for a firing that
a correctly attached recorder observed, because `enable!` always precedes
`fire!`.

# Fields
 - `clock::K` — the key of the clock that fired.
 - `when::T` — absolute firing time.
 - `te::T` — the enabling time, the zero-point of `distribution`; `te ≤ when`.
 - `distribution::Distributions.UnivariateDistribution` — the lifetime
   distribution current at firing, retained so the record replays offline
   without re-deriving which distribution was live.
 - `u::Float64` — the total-lifetime survival uniform, `ccdf(distribution, when - te)`.
 - `logu::Float64` — `log(u)`, the log-survival coordinate the identity inverts.
"""
struct ClockFiredRecord{K,T}
    clock::K
    when::T
    te::T
    distribution::UnivariateDistribution
    u::Float64
    logu::Float64
end


"""
    TrajectoryRecorder{K,T}()

A watcher that records the ordered sequence of clock firings of a trajectory.

It maintains the live `(distribution, te)` table exactly as [`TrackWatcher`](@ref)
does — inheriting `enable!`/`disable!` from `CompetingClocks.EnabledWatcher`
— and overrides `fire!` to append a [`ClockFiredRecord`](@ref). At each firing it
back-calculates the survival-space uniform in log space,

```julia
logu = logccdf(distribution, when - te)   # ⇒  when == te + invlogccdf(distribution, logu)
u    = exp(logu)
```

so the recorder never asks the sampler for randomness: it is **sampler-agnostic**
and produces the same record on `FirstReaction`, `CombinedNextReaction`, or any
other correct sampler. The recorded identity is a CONTRACT OBLIGATION every
firing must satisfy, which is what makes the record a sufficient statistic for
offline likelihood replay and pathwise (IPA) derivatives.

Attach one to a running `SamplingContext` with [`with_recorder`](@ref), or drive
it directly like any watcher. Read the firings back with [`recorded_firings`](@ref)
and stamp the observation horizon (the time censoring terms extend to) with
[`close_record!`](@ref).

```julia
rec = TrajectoryRecorder{Int,Float64}()
enable!(rec, 3, Weibull(1.5, 2.0), 0.0, 0.0)
fire!(rec, 3, 1.7)              # records (3, 1.7, 0.0, dist, u, logu)
close_record!(rec, 10.0)       # observation ran to horizon 10.0
recorded_firings(rec)          # Vector{ClockFiredRecord{Int,Float64}}
```
"""
mutable struct TrajectoryRecorder{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    firings::Vector{ClockFiredRecord{K,T}}
    horizon::T
    closed::Bool
    TrajectoryRecorder{K,T}() where {K,T} =
        new(Dict{K,EnablingEntry{K,T}}(), ClockFiredRecord{K,T}[], zero(T), false)
end


# A clone is an empty recorder of the same type, matching TrackWatcher and
# DebugWatcher: `clone` seeds a fresh run, it does not duplicate a history.
clone(rec::TrajectoryRecorder{K,T}) where {K,T} = TrajectoryRecorder{K,T}()


# reset! must clear the firing history as well as the live table; the generic
# EnabledWatcher reset! only empties `enabled`, which would leak firings from a
# previous run into the next.
function reset!(rec::TrajectoryRecorder{K,T}) where {K,T}
    empty!(rec.enabled)
    empty!(rec.firings)
    rec.horizon = zero(T)
    rec.closed = false
    return nothing
end


# copy_clocks! carries BOTH the live table and the firing history so a split /
# branch continuation inherits the path recorded so far, not just the enabled
# set. The generic EnabledWatcher copy_clocks! copies only `enabled`.
function copy_clocks!(dst::TrajectoryRecorder{K,T}, src::TrajectoryRecorder{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    copy!(dst.firings, src.firings)
    dst.horizon = src.horizon
    dst.closed = src.closed
    return dst
end


"""
    fire!(recorder::TrajectoryRecorder, clock, when)

Record that `clock` fired at time `when`. Looks up the clock's live
`(distribution, te)`, back-calculates the survival-space uniform in log space,
appends a [`ClockFiredRecord`](@ref), and removes the clock from the live table.

The stored `u` is the uniform of the clock's TOTAL lifetime,
`ccdf(distribution, when - te)`, NOT a uniform conditioned on survival to the
enabling call. For a left-shifted enabling (`te < when`) this distinction
matters: the identity is anchored at `te`, so a clock aged before it was
observed still records the total-lifetime `u`.
"""
function fire!(rec::TrajectoryRecorder{K,T}, clock::K, when::T) where {K,T}
    entry = get(rec.enabled, clock, nothing)
    if entry === nothing
        # enable! always precedes fire! for an attached recorder, so this branch
        # is a defensive degradation, not a supported path: no (dist, te) means
        # no back-calculation, so u/logu are NaN and the distribution is a
        # never-firing placeholder.
        push!(rec.firings, ClockFiredRecord{K,T}(clock, when, when, Never(), NaN, NaN))
        return nothing
    end
    dist = entry.distribution
    te = entry.te
    # The retained-draw identity lives in log-survival space:
    #   logu = logccdf(dist, when - te)  ⇒  when == te + invlogccdf(dist, logu).
    # Deriving logu first (rather than u = ccdf(...)) keeps a deep tail from
    # underflowing u to 0.0, which would make the inversion lose `when`.
    logu = Float64(logccdf(dist, when - te))
    u = exp(logu)
    push!(rec.firings, ClockFiredRecord{K,T}(clock, when, te, dist, u, logu))
    delete!(rec.enabled, clock)
    return nothing
end


"""
    close_record!(recorder, horizon)

Stamp the observation end. `horizon` is the time the trajectory was observed to,
and it is where the survival (censoring) terms of every clock still enabled at
the end of the run extend to. Marks the recorder closed and returns it.
"""
function close_record!(rec::TrajectoryRecorder{K,T}, horizon) where {K,T}
    rec.horizon = convert(T, horizon)
    rec.closed = true
    return rec
end


"""
    recorded_firings(recorder)

The ordered `Vector{ClockFiredRecord{K,T}}` of firings observed so far.
"""
recorded_firings(rec::TrajectoryRecorder) = rec.firings


"""
    horizon(recorder)

The observation horizon stamped by [`close_record!`](@ref), or `zero(T)` if the
record has not been closed. Check [`isclosed`](@ref) to distinguish an unclosed
record from one closed at time zero.
"""
horizon(rec::TrajectoryRecorder) = rec.horizon

"""
    isclosed(recorder)

Whether [`close_record!`](@ref) has stamped an observation horizon.
"""
isclosed(rec::TrajectoryRecorder) = rec.closed


# ---------------------------------------------------------------------------
# Attaching a recorder to a running SamplingContext.
#
# The context fans every verb out to a compile-time tuple of watchers. There is
# no builder flag for a bespoke watcher, so a recorder joins by rebuilding the
# context around an extended watcher tuple. The rebuilt context SHARES the
# sampler, rng, and delayed state with the original (it is a live view with
# one more observer), so use only the returned context afterward.
#
# The recorder is intentionally NOT a `LikelihoodWatcher`: on a vector-enable it
# must see the SAMPLED distribution (the one the sampler actually drew from),
# not the whole importance-sampling vector, because `u` is back-calculated
# against the distribution that produced the firing. The default fan-out already
# hands it the sampled distribution, so no change to `context.jl` is needed.
# ---------------------------------------------------------------------------

"""
    attach_watcher(ctx::SamplingContext, w) -> SamplingContext

Return a `SamplingContext` identical to `ctx` but with watcher `w` appended to
its watcher tuple. The result shares `ctx`'s sampler, rng, and delayed
state, so drive only the returned context after calling this.
"""
function attach_watcher(ctx::SamplingContext{K,T,Sampler,RNG,W,DS}, w) where {K,T,Sampler,RNG,W,DS}
    watchers = (ctx.watchers..., w)
    return SamplingContext{K,T,Sampler,RNG,typeof(watchers),DS}(
        ctx.sampler, ctx.rng, watchers,
        ctx.split_weight, ctx.time, ctx.fixed_start, ctx.sample_distribution, ctx.delayed,
    )
end


"""
    with_recorder(ctx::SamplingContext) -> (ctx′, recorder)

Build a [`TrajectoryRecorder`](@ref) keyed to the context's internal clock/time
types, attach it to `ctx`, and return the recorder together with the extended
context `ctx′`. Drive `ctx′` from then on; read firings from `recorder`.

```julia
ctx = SamplingContext(SamplerBuilder(K, T; method=NextReactionMethod()), rng)
ctx, rec = with_recorder(ctx)
# ... run the simulation on ctx ...
close_record!(rec, horizon)
recorded_firings(rec)
```
"""
function with_recorder(ctx::SamplingContext)
    Kint = keytype(ctx.sampler)
    Tt = timetype(ctx.sampler)
    rec = TrajectoryRecorder{Kint,Tt}()
    return attach_watcher(ctx, rec), rec
end
