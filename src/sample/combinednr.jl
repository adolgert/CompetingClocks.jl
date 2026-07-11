
using DataStructures: MutableBinaryHeap, update!

"""
This function decides whether a particular distribution can be sampled faster
and more accurately using its cumulative distribution function or using
the log of its cumulative distribution function, also called the integrated
hazard. The former is used for the Next Reaction method by Gibson and Bruck.
The latter is used by the Modified Next Reaction method of Anderson.
We are calling the first a linear space and the second a logarithmic space.
"""
function sampling_space(x::T) where {T<:UnivariateDistribution}
    sampling_space(typeof(x))::Union{Type{LinearSampling},Type{LogSampling}}
end
abstract type SamplingSpaceType end
struct LinearSampling <: SamplingSpaceType end
struct LogSampling <: SamplingSpaceType end
sampling_space(::Type) = LinearSampling
sampling_space(::Type{<:Distributions.Arcsine}) = LinearSampling
sampling_space(::Type{<:Distributions.BetaPrime}) = LinearSampling
sampling_space(::Type{<:Distributions.Biweight}) = LinearSampling
sampling_space(::Type{<:Distributions.Beta}) = LinearSampling
sampling_space(::Type{<:Distributions.Cauchy}) = LinearSampling
sampling_space(::Type{<:Distributions.Cosine}) = LinearSampling
sampling_space(::Type{<:Distributions.Epanechnikov}) = LinearSampling
sampling_space(::Type{<:Distributions.Erlang}) = LogSampling
sampling_space(::Type{<:Distributions.Exponential}) = LogSampling
sampling_space(::Type{<:Distributions.Frechet}) = LinearSampling
sampling_space(::Type{<:Distributions.Gamma}) = LogSampling
sampling_space(::Type{<:Distributions.GeneralizedPareto}) = LinearSampling
sampling_space(::Type{<:Distributions.Gumbel}) = LinearSampling
sampling_space(::Type{<:Distributions.InverseGamma}) = LinearSampling
sampling_space(::Type{<:Distributions.InverseGaussian}) = LinearSampling
sampling_space(::Type{<:Distributions.JohnsonSU}) = LinearSampling
sampling_space(::Type{<:Distributions.Kolmogorov}) = LinearSampling
sampling_space(::Type{<:Distributions.Kumaraswamy}) = LinearSampling
sampling_space(::Type{<:Distributions.Laplace}) = LogSampling
sampling_space(::Type{<:Distributions.Levy}) = LinearSampling
sampling_space(::Type{<:Distributions.Lindley}) = LinearSampling
sampling_space(::Type{<:Distributions.Logistic}) = LinearSampling
sampling_space(::Type{<:Distributions.LogitNormal}) = LinearSampling
sampling_space(::Type{<:Distributions.LogNormal}) = LinearSampling
sampling_space(::Type{<:Distributions.Normal}) = LinearSampling
sampling_space(::Type{<:Distributions.NormalCanon}) = LinearSampling
sampling_space(::Type{<:Distributions.Pareto}) = LinearSampling
sampling_space(::Type{<:Distributions.PGeneralizedGaussian}) = LinearSampling
sampling_space(::Type{<:Distributions.Rayleigh}) = LinearSampling
sampling_space(::Type{<:Distributions.Rician}) = LinearSampling
sampling_space(::Type{<:Distributions.SkewedExponentialPower}) = LinearSampling
sampling_space(::Type{<:Distributions.SymTriangularDist}) = LinearSampling
sampling_space(::Type{<:Distributions.Triweight}) = LinearSampling
sampling_space(::Type{<:Distributions.Uniform}) = LinearSampling
sampling_space(::Type{<:Distributions.Weibull}) = LogSampling

# The following four support functions are used by the CombinedNextReaction
# sampler, and their use is decided by the `sampling_space()` above.
# I had some trouble getting these functions to be type stable, and their lack
# of type stability hurt performance. Test performance by using
# tests/time_combinednr.jl.
function get_survival_zero(::T) where {T<:UnivariateDistribution}
    get_survival_zero(sampling_space(T))::Float64
end
function get_survival_zero(::Type{T}) where {T<:UnivariateDistribution}
    get_survival_zero(sampling_space(T))::Float64
end
get_survival_zero(::Type{LinearSampling}) = 0.0::Float64
get_survival_zero(::Type{LogSampling}) = -Inf::Float64

draw_space(::Type{LinearSampling}, rng) = rand(rng, Uniform())
draw_space(::Type{LogSampling}, rng) = rand(rng, Exponential())

function survival_space(::Type{T}, dist, sample) where {T<:UnivariateDistribution}
    survival_space(sampling_space(T), dist, sample)
end
survival_space(::Type{LinearSampling}, dist, sample) = ccdf(dist, sample)::Float64
survival_space(::Type{LogSampling}, dist, sample) = logccdf(dist, sample)::Float64

function invert_space(::Type{T}, dist, survival) where {T<:UnivariateDistribution}
    invert_space(sampling_space(T), dist, survival)
end
invert_space(::Type{LinearSampling}, dist, survival) = cquantile(dist, survival)::Float64
invert_space(::Type{LogSampling}, dist, survival) = invlogccdf(dist, survival)::Float64


struct NRTransition{T}
    heap_handle::Int
    survival::T # value of S_j or Λ_j
    distribution::UnivariateDistribution
    te::T  # Enabling time of distribution
    t0::T  # Enabling time of transition
end


"""
    CombinedNextReaction{KeyType,TimeType}(seed; coupling=:carry)

This combines Next Reaction Method and Modified Next Reaction Method.
The Next Reaction Method is from Gibson and Bruck in their 2000 paper called
["Efficient Exact Stochastic Simulation of Chemical Systems with Many Species
and Many Channels"](https://doi.org/10.1021/jp993732q). 
The Modified Next Reaction Method is from David F. Anderson's 2007 paper, 
["A modified Next Reaction Method for simulating chemical systems with time dependent propensities and delays"](https://doi.org/10.1063/1.2799998). 
Both methods reuse draws of random numbers. The former works by accumulating 
survival of a distribution in a linear space and the latter works by accumulating 
survival of a distribution in a log space.

Each enabled clock specifies a univariate distribution from the `Distributions`
package. Every distribution is more precise being sampled in the manner
of the Next Reaction method (linear space) or the manner of the Modified
Next Reaction method (log space). This sampler
chooses which space to use depending on the
type of the `UnivariateDistribution` and based on performance timings that
are done during package testing. Defaults are set for those distributions
included in the `Distributions.jl` package. If you want to add a distribution,
then define:

```julia
sampling_space(::MyDistribution) = LogSampling
```

If you want to override a choice in the library, then create a sub-type of the
given distribution, and specify its sampling space.

```julia
struct LinearGamma <: Distributions.Gamma end
sampling_space(::LinearGamma) = LinearSampling
```

If you want to test a distribution, look at `tests/nrmetric.jl` to see how
distributions are timed.

The `coupling` keyword (default `:carry`) fixes, at construction, how
[`reenable!`](@ref) realizes a mid-flight distribution change: `:carry` maps
the retained survival draw through the change deterministically (the
Anderson–Kurtz random-time-change coupling; consumes no randomness and is what
pathwise/IPA derivatives need), while `:redraw` draws the remaining lifetime
fresh, conditioned on age. Both agree in law. Read it back with
[`coupling`](@ref). The `:carry` default restores this sampler's historical
behavior on re-enabling an enabled clock.
"""
mutable struct CombinedNextReaction{K,T} <: SSA{K,T}
    firing_queue::MutableBinaryHeap{OrderedSample{K,T},Base.ForwardOrdering}
    transition_entry::Dict{K,NRTransition{T}}
    streams::KeyedStreams{K}
    coupling::Symbol
end


function CombinedNextReaction{K,T}(seed=_DEFAULT_STREAM_SEED; coupling::Symbol=:carry) where {K,T<:ContinuousTime}
    validate_coupling(CombinedNextReaction, coupling)
    heap = MutableBinaryHeap{OrderedSample{K,T},Base.ForwardOrdering}()
    CombinedNextReaction{K,T}(heap, Dict{K,NRTransition{T}}(), KeyedStreams{K}(seed), coupling)
end

coupling(nr::CombinedNextReaction) = nr.coupling

# A full-state clone copies the firing queue, the retained-survival table, AND
# the keyed streams (generator states, counts). Because this sampler REUSES a
# clock's survival draw across enabling episodes, carrying the stream states is
# what keeps the clone's future draws identical to the original's — the coupling
# the CRN successor relies on.
function clone(nr::CombinedNextReaction{K,T}) where {K,T}
    c = CombinedNextReaction{K,T}(nr.streams.seed; coupling=nr.coupling)
    c.firing_queue = deepcopy(nr.firing_queue)
    copy!(c.transition_entry, nr.transition_entry)
    c.streams = copy(nr.streams)
    return c
end

similar_sampler(nr::CombinedNextReaction{K,T}) where {K,T} =
    CombinedNextReaction{K,T}(nr.streams.seed; coupling=nr.coupling)

rekey_streams!(nr::CombinedNextReaction, seed) = (rekey_streams!(nr.streams, seed); nr)

function reset!(nr::CombinedNextReaction)
    empty!(nr.firing_queue)
    @assert isempty(nr.firing_queue)
    empty!(nr.transition_entry)
    nothing
end

function copy_clocks!(dst::CombinedNextReaction{K,T}, src::CombinedNextReaction{K,T}) where {K,T}
    dst.firing_queue = deepcopy(src.firing_queue)
    copy!(dst.transition_entry, src.transition_entry)
    dst.streams = copy(src.streams)
    # The coupling is a construction-time property, but copy_clocks! promises a
    # full replacement of the destination's state, so the destination must
    # re-evaluate clocks the same way the source would have.
    dst.coupling = src.coupling
    return dst
end


function _jitter_space(nr::CombinedNextReaction{K,T}, clock, record, ::Type{S}, when::T,
        ) where {K,T,S<:SamplingSpaceType}
    if record.heap_handle > 0
        tau, shift_survival = sample_shifted(stream_for!(nr.streams, clock), record.distribution, S, record.te, when)
        sample = OrderedSample{K,T}(clock, tau)
        update!(nr.firing_queue, record.heap_handle, sample)
        nr.transition_entry[clock] = NRTransition{T}(
            record.heap_handle, shift_survival, record.distribution, record.te, when
        )
    else
        # A retained-disabled entry (heap_handle == 0) banks residual survival
        # for Anderson/Gibson-Bruck reuse; decorrelation must extinguish that
        # bank too, or a re-enabled clock replays pre-jitter randomness.
        # NRTransition is immutable, so REPLACE the entry — assigning through
        # `record.survival = ...` threw whenever a disabled entry existed,
        # which made jitter! unusable on any sampler that had ever disabled a
        # clock (found by ClockGradients' branchable-world rekey sweep).
        nr.transition_entry[clock] = NRTransition{T}(
            0, get_survival_zero(S), record.distribution, record.te, record.t0
        )
    end
end


function jitter!(nr::CombinedNextReaction{K,T}, when::T) where {K,T}
    for (clock, record) in nr.transition_entry
        _jitter_space(nr, clock, record, sampling_space(record.distribution), when)
    end
end


@doc raw"""
For the first reaction sampler, you can call next() multiple times and get
different, valid, answers. For a CombinedNextReaction sampler, next() is
non-mutating: it reads the minimum entry of the firing queue and returns it
without changing any internal state and without consuming any random numbers.
Repeated calls to next(), with no intervening enable!, disable!, or fire!,
return the same cached reservation. Committing to a firing is done separately
by calling fire!.
"""
function next(nr::CombinedNextReaction{K,T}, when::T) where {K,T}
    if !isempty(nr.firing_queue)
        least = first(nr.firing_queue)
        return (least.time, least.key)
    else
        # Return type is Tuple{Float64, Union{Nothing,T}} because T is not default-constructible.
        return (typemax(T), nothing)
    end
end


# Implementation note: This function and others below are parametrized on the
# sampling space. If you parametrize on individual distributions, it will create
# too many specializations, so using the SamplingSpaceType is a nice compromise.
function sample_shifted(
    rng::AbstractRNG,
    distribution::UnivariateDistribution,
    ::Type{S},
    te::T,
    when::T
) where {S<:SamplingSpaceType,T<:ContinuousTime}
    if te < when
        shifted_distribution = truncated(distribution, when - te, typemax(T))
        sample = rand(rng, shifted_distribution)
        tau = te + sample
        survival = survival_space(S, shifted_distribution, sample)
    else  # te >= when
        # The distribution starts in the future
        sample = rand(rng, distribution)
        tau = te + sample
        survival = survival_space(S, distribution, sample)
    end
    (tau, survival)
end


function sample_by_inversion(
    distribution::UnivariateDistribution, ::Type{S}, te::T, when::T, survival::T
) where {S<:SamplingSpaceType,T<:ContinuousTime}
    if te < when
        te + invert_space(S, truncated(distribution, when - te, typemax(T)), survival)
    else   # te > when
        te + invert_space(S, distribution, survival)
    end
end


@doc raw"""
This updates the survival for a transition in the linear space, according to
Gibson and Bruck.
Transition was enabled between time record t_0 and t_n.
Divide the survival by the conditional survival between t_0 and t_n.
t_e can be before t_0, at t_0, between t_0 and t_n, or at t_n, or after t_n.

``u=\exp\left(-\int_{t_e}^{t_n}\lambda_0(s-t_e)ds\right)\exp\left(-\int_{t_n}^{\tau'}\lambda_{n}(s-t_e)ds\right)``

"""
function consume_survival(
    record::NRTransition, distribution::UnivariateDistribution, ::Type{S}, tn::T
) where {S<:LinearSampling,T<:ContinuousTime}
    survive_te_tn = if record.te < tn
        ccdf(distribution, tn - record.te)::T
    else
        one(T)
    end
    survive_te_t0 = if record.te < record.t0
        ccdf(distribution, record.t0 - record.te)::T
    else
        one(T)
    end
    record.survival * (survive_te_t0 / survive_te_tn)
end


@doc raw"""
This updates the survival for a transition in log space, according to
Anderson's method.

``\ln u=-\int_{t_e}^{t_n}\lambda_0(s-t_e)ds - \int_{t_n}^{\tau'}\lambda_{n}(s-t_e)ds``

"""
function consume_survival(
    record::NRTransition, distribution::UnivariateDistribution, ::Type{S}, tn::T
) where {S<:LogSampling,T<:ContinuousTime}
    log_survive_te_tn = if record.te < tn
        logccdf(distribution, tn - record.te)::T
    else
        zero(T)
    end
    log_survive_te_t0 = if record.te < record.t0
        logccdf(distribution, record.t0 - record.te)::T
    else
        zero(T)
    end
    record.survival - (log_survive_te_tn - log_survive_te_t0)
end


function enable!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T) where {K,T}
    enable!(nr, clock, distribution, sampling_space(distribution), te, when)
    nothing
end


function enable!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution, ::Type{S},
    te::T, when::T) where {K,T,S<:SamplingSpaceType}

    # Three cases: a) never been enabled b) currently enabled c) was disabled.
    record = get(
        nr.transition_entry,
        clock,
        NRTransition{T}(0, get_survival_zero(S), Never(), zero(T), zero(T))
    )
    heap_handle = record.heap_handle

    # if the transition needs to be re-drawn.
    if record.survival <= get_survival_zero(S)
        tau, shift_survival = sample_shifted(stream_for!(nr.streams, clock), distribution, S, te, when)
        sample = OrderedSample{K,T}(clock, tau)
        if record.heap_handle > 0
            update!(nr.firing_queue, record.heap_handle, sample)
        else
            heap_handle = push!(nr.firing_queue, sample)
        end
        nr.transition_entry[clock] = NRTransition{T}(
            heap_handle, shift_survival, distribution, te, when
        )

        # The transition has remaining lifetime.
    else
        # The transition was previously enabled.
        if record.heap_handle > 0
            # Consider te the same if the mantissa is within 2 bits of precision.
            same_te = abs(te - record.te) < 2 * eps(te)
            if same_te && distribution == record.distribution
                # No change. It's common to re-enable an already-enabled distribution.
            else
                # Account for time between when this was last enabled and now.
                survival_remain = consume_survival(record, record.distribution, S, when)
                tau = sample_by_inversion(distribution, S, te, when, survival_remain)
                entry = OrderedSample{K,T}(clock, tau)
                update!(nr.firing_queue, record.heap_handle, entry)
                nr.transition_entry[clock] = NRTransition{T}(
                    heap_handle, survival_remain, distribution, te, when
                )
            end

            # The transition was previously disabled.
        else
            tau = sample_by_inversion(distribution, S, te, when, record.survival)
            heap_handle = push!(nr.firing_queue, OrderedSample{K,T}(clock, tau))
            nr.transition_entry[clock] = NRTransition{T}(
                heap_handle, record.survival, distribution, te, when
            )
        end
    end
    nothing
end


@doc raw"""
    fire!(nr::CombinedNextReaction, clock, when)

Commit to firing `clock` at time `when`. This removes the clock from the firing
queue and consumes its draw completely by setting its stored survival to the
zero of its sampling space (`0.0` for linear space, `-Inf` for log space). As a
result, re-enabling this clock later takes the fresh-redraw branch of `enable!`,
which is the correct behavior for a draw that has been fully consumed by firing.

This differs from `disable!`, which preserves the clock's *remaining* survival
so that a later re-enable can reuse the draw (the Anderson/Gibson-Bruck
draw-reuse property).
"""
function fire!(nr::CombinedNextReaction{K,T}, clock::K, when::T) where {K,T<:ContinuousTime}
    record = nr.transition_entry[clock]
    delete!(nr.firing_queue, record.heap_handle)
    nr.transition_entry[clock] = NRTransition{T}(
        0,
        get_survival_zero(sampling_space(record.distribution)),
        record.distribution,
        record.te,
        when
    )
    nothing
end


function disable!(nr::CombinedNextReaction{K,T}, clock::K, when::T) where {K,T<:ContinuousTime}
    record = nr.transition_entry[clock]
    delete!(nr.firing_queue, record.heap_handle)
    nr.transition_entry[clock] = NRTransition{T}(
        0,
        consume_survival(record, record.distribution, sampling_space(record.distribution), when),
        record.distribution,
        record.te,
        when
    )
    nothing
end


function steploglikelihood(nr::CombinedNextReaction, t0, t, which_fires)
    # We need to adapt our dictionary into a named tuple for consumption by the
    # calculator of the step log-likelihood. Fired and disabled clocks are retained
    # in transition_entry (with heap_handle == 0) for draw reuse on re-enable, so we
    # filter to enabled clocks (heap_handle > 0) to exclude their spurious survival.
    return _steploglikelihood(
        ((clock=k, distribution=v.distribution, te=v.te)
         for (k, v) in pairs(nr.transition_entry) if v.heap_handle > 0),
        t0,
        t,
        which_fires
    )
end


"""
    getindex(sampler::CombinedNextReaction{K,T}, clock::K)

For the `CombinedNextReaction` sampler, returns the stored firing time associated to the clock.
"""
function Base.getindex(nr::CombinedNextReaction{K,T}, clock::K) where {K,T}
    if haskey(nr.transition_entry, clock)
        heap_handle = getfield(nr.transition_entry[clock], :heap_handle)
        return getfield(nr.firing_queue[heap_handle], :time)
    else
        throw(KeyError(clock))
    end
end

function Base.keys(nr::CombinedNextReaction)
    return collect(keys(nr.transition_entry))
end

function Base.length(nr::CombinedNextReaction)
    return length(nr.transition_entry)
end


function isenabled(nr::CombinedNextReaction, clock)
    haskey(nr.transition_entry, clock) && nr.transition_entry[clock].heap_handle > 0
end


# A set of all enabled clock keys for a CombinedNextReaction method.
# We make a custom Set implementation because the information is in the
# CombinedNextReaction object, but it's spread across a Heap and a Dictionary.
# This helper class should make it much more efficient to iterate the set.
struct NextReactionEnabled{C,T,K} <: AbstractSet{C}
    nr::T
    keys::K
end

_has_handle(nre::NextReactionEnabled, key) = nre.nr.transition_entry[key].heap_handle > 0

function Base.iterate(nre::NextReactionEnabled)
    res = iterate(nre.keys)
    res === nothing && return res
    while !_has_handle(nre, res[1])
        res = iterate(nre.keys, res[2])
        res === nothing && return res
    end
    return res
end


function Base.iterate(nre::NextReactionEnabled, state)
    res = iterate(nre.keys, state)
    res === nothing && return res
    while !_has_handle(nre, res[1])
        res = iterate(nre.keys, res[2])
        res === nothing && return res
    end
    return res
end

Base.length(nre::NextReactionEnabled) = length(nre.nr.firing_queue)
Base.in(x, nre::NextReactionEnabled) = isenabled(nre.nr, x)
Base.eltype(::Type{NextReactionEnabled{C}}) where {C} = C


function enabled(nr::CombinedNextReaction{K,T}) where {K,T}
    kks = keys(nr.transition_entry)
    NextReactionEnabled{K,typeof(nr),typeof(kks)}(nr, kks)
end


function Base.haskey(nr::CombinedNextReaction{K,T}, clock::K) where {K,T}
    return haskey(nr.transition_entry, clock) && nr.transition_entry[clock].heap_handle > 0
end

Base.haskey(dc::CombinedNextReaction{K,T}, clock) where {K,T} = false


# --- estimator-facing verbs -------------------------------------------------
# CombinedNextReaction is the one sampler in this package that retains draw
# randomness across state changes (the Anderson/Gibson-Bruck reuse property), so
# it is the only one that supports `retained_draw`. It is a scheduling backend
# with a per-clock enabling-time table, so it also supports `enabled_ages` and
# `force_fire!`. Disabled and fired clocks are kept in `transition_entry` with
# `heap_handle == 0` for draw reuse; every estimator verb here excludes them.

supports_enabled_ages(::Type{<:CombinedNextReaction}) = true
supports_force(::Type{<:CombinedNextReaction}) = true
supports_retained_draw(::Type{<:CombinedNextReaction}) = true
# CombinedNextReaction retains each clock's survival across state changes, so it
# can carry a mid-flight distribution change deterministically — indeed its
# historical enable!-on-enabled-key behavior ALREADY IS carry (see reenable!).
supports_carry(::Type{<:CombinedNextReaction}) = true


"""
    reenable!(nr::CombinedNextReaction, clock, distribution, te, when)

Re-evaluate `clock`'s distribution mid-flight. Which pathwise coupling realizes
the change is the sampler's construction-time [`coupling`](@ref) field.

`:carry` reuses the retained survival: this is EXACTLY the already-enabled branch
of [`enable!`](@ref) (`consume_survival` re-references the stored survival from
its last reschedule to `when`, then `sample_by_inversion` maps it into the new
distribution), so CombinedNextReaction's historical enable!-on-rate-change
behavior has always been carry. Carry consumes no randomness and, for an
identical distribution, leaves the schedule bit-for-bit unchanged.

`:redraw` discards the retained survival and draws the remaining lifetime fresh
from the clock's OWN keyed stream, conditioned on the current age (via
`sample_shifted`, which truncates at `when - te`). This is the same fresh draw a
plain first-enable would make, so under `:redraw` even an identical distribution
generally moves the schedule.
"""
function reenable!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T) where {K,T}
    haskey(nr, clock) || throw(ArgumentError(
        "reenable! needs clock $clock currently enabled (heap_handle > 0); " *
        "use enable! to start a clock that is not enabled."))
    # The constructor validated the field, so only the two couplings reach here.
    if nr.coupling === :carry
        # The already-enabled branch of enable! is the carry map.
        enable!(nr, clock, distribution, te, when)
    else
        _reenable_redraw!(nr, clock, distribution, sampling_space(distribution), te, when)
    end
    nothing
end


# Redraw coupling: the fresh-draw branch of enable! (sample_shifted from te at
# the current time), taken unconditionally even when the distribution is
# unchanged. The clock keeps its heap slot; only its scheduled time and stored
# survival change, drawn from the clock's own keyed stream so a coupled clone
# pair redraws identically.
function _reenable_redraw!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    ::Type{S}, te::T, when::T) where {K,T,S<:SamplingSpaceType}
    record = nr.transition_entry[clock]
    tau, shift_survival = sample_shifted(stream_for!(nr.streams, clock), distribution, S, te, when)
    update!(nr.firing_queue, record.heap_handle, OrderedSample{K,T}(clock, tau))
    nr.transition_entry[clock] = NRTransition{T}(
        record.heap_handle, shift_survival, distribution, te, when
    )
    nothing
end

# Only currently-enabled clocks (heap_handle > 0). A disabled or fired entry is
# retained for draw reuse but is not part of the enabled set, so it must not
# appear among the ages an estimator races.
enabling_times(nr::CombinedNextReaction) =
    ((clock, rec.te) for (clock, rec) in nr.transition_entry if rec.heap_handle > 0)


# The log-survival coordinate of a stored survival value. CombinedNextReaction
# keeps survival in one of two sampling spaces per distribution family (see
# `sampling_space`): linear space stores the survival probability S_j itself,
# log space stores the integrated hazard Λ_j = log S_j. Normalizing both to the
# log coordinate is what lets `retained_draw` report one identity for every
# family.
_log_survival(::Type{LinearSampling}, survival) = log(survival)
_log_survival(::Type{LogSampling}, survival) = survival


"""
    retained_draw(nr::CombinedNextReaction, clock) -> (te=te, u=u, logu=logu)

Return the enabling time and the survival-space uniform behind `clock`'s current
tentative firing time, satisfying `tentative == te + invlogccdf(dist, logu)`.

CombinedNextReaction stores its per-clock survival in whichever sampling space
(`LinearSampling` or `LogSampling`) suits the distribution family, and, for a
left-shifted enabling (`te < t0`, where `t0` is the last (re)schedule time),
that survival is measured within the distribution TRUNCATED at the age
`t0 - te`. This normalizes both to the log-survival coordinate of the clock's
TOTAL lifetime from `te`:

    logu = _log_survival(space, survival) + (te < t0 ? logccdf(dist, t0 - te) : 0)

so the same identity holds whatever the family or the enabling shift.
"""
function retained_draw(nr::CombinedNextReaction{K,T}, clock::K) where {K,T}
    haskey(nr.transition_entry, clock) || throw(KeyError(clock))
    rec = nr.transition_entry[clock]
    rec.heap_handle > 0 || throw(ArgumentError(
        "retained_draw needs clock $clock currently enabled (heap_handle > 0); " *
        "a disabled or fired clock has no tentative firing time."))
    S = sampling_space(rec.distribution)
    shift = rec.te < rec.t0 ? logccdf(rec.distribution, rec.t0 - rec.te) : zero(T)
    logu = _log_survival(S, rec.survival) + shift
    return (te=rec.te, u=exp(logu), logu=logu)
end


"""
    force_fire!(nr::CombinedNextReaction, clock, tstar)

Fire `clock` at `tstar` regardless of the race. The fired clock is consumed
exactly as `fire!` consumes it — its survival is zeroed so a later re-enable
redraws fresh, NOT preserved as `disable!` would. Each surviving enabled clock
is repaired with keep-if-later / redraw-if-passed: a stored schedule already
past `tstar` is kept (conditioned on the race being decided at `tstar`, the
unconditional draw that exceeds `tstar` is already survival-conditioned), while
a schedule `≤ tstar` — a promise the force overran — is redrawn from the
lifetime truncated past `tstar` via `sample_shifted`, drawing from that loser's
OWN keyed stream so the redraw is identical across a coupled clone pair, and its
stored survival and schedule updated so the retained-draw identity stays true.

See the unbiasedness precondition on the generic `force_fire!`.
"""
function force_fire!(nr::CombinedNextReaction{K,T}, clock::K, tstar::T) where {K,T<:ContinuousTime}
    haskey(nr, clock) || throw(KeyError(clock))  # haskey requires heap_handle > 0
    fired = nr.transition_entry[clock]
    Sfired = sampling_space(fired.distribution)
    delete!(nr.firing_queue, fired.heap_handle)
    # Consume the draw like fire!, not disable!: zero the survival so re-enabling
    # takes the fresh-redraw branch.
    nr.transition_entry[clock] = NRTransition{T}(
        0, get_survival_zero(Sfired), fired.distribution, fired.te, tstar
    )
    # Repair survivors. Collect keys first so reassigning dict entries in the
    # loop cannot disturb iteration.
    for k in collect(keys(nr.transition_entry))
        rec = nr.transition_entry[k]
        rec.heap_handle > 0 || continue  # skip disabled/fired (incl. the just-fired clock)
        scheduled = getfield(nr.firing_queue[rec.heap_handle], :time)
        scheduled > tstar && continue    # keep-if-later
        # Redraw-if-passed: fresh draw from the age-conditioned law past tstar.
        # Storing t0 = tstar makes tstar the truncation reference, so the redrawn
        # survival reconstructs through retained_draw exactly as an enable would.
        S = sampling_space(rec.distribution)
        tau, shift_survival = sample_shifted(stream_for!(nr.streams, k), rec.distribution, S, rec.te, tstar)
        update!(nr.firing_queue, rec.heap_handle, OrderedSample{K,T}(k, tau))
        nr.transition_entry[k] = NRTransition{T}(
            rec.heap_handle, shift_survival, rec.distribution, rec.te, tstar
        )
    end
    nothing
end
