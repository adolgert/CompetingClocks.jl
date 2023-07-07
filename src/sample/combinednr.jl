
using DataStructures

export sampling_space
export CombinedNextReaction

"""
This function decides whether a particular distribution can be sampled faster
and more accurately using its cumulative distribution function or using
the log of its cumulative distribution function, also called the integrated
hazard. The former is used for the Next Reaction method by Gibson and Bruck.
The latter is used by the Modified Next Reaction method of Anderson.
We are calling the first a linear space and the second a logarithmic space.
"""
sampling_space(x) = sampling_space(typeof(x))
abstract type SamplingSpaceType end
struct LinearSampling <: SamplingSpaceType end
struct LogSampling <: SamplingSpaceType end
sampling_space(::Type) = LinearSampling
sampling_space(::Type{Distributions.Exponential}) = LogSampling
sampling_space(::Type{Distributions.Frechet}) = LogSampling
sampling_space(::Type{Distributions.Gamma}) = LogSampling
sampling_space(::Type{Distributions.Semicircle}) = LinearSampling

# The following four support functions are used by the CombinedNextReaction
# sampler, and their use is decided by the `sampling_space()` above.
get_survival_zero(x::T) where {T} = get_survival_zero(sampling_space(T))
get_survival_zero(::Type{LinearSampling}) = 0.0
get_survival_zero(::Type{LogSampling}) = -Inf

get_survival_max(x::T) where {T} = get_survival_max(sampling_space(T))
get_survival_max(::Type{LinearSampling}) = 1.0
get_survival_max(::Type{LogSampling}) = 0.0

survival_space(dist, sample) = survival_space(sampling_space(dist), dist, sample)
survival_space(::Type{LinearSampling}, dist, sample) = ccdf(dist, sample)
survival_space(::Type{LogSampling}, dist, sample) = logccdf(dist, sample)

invert_space(dist, survival) = invert_space(sampling_space(dist), dist, survival)
invert_space(::Type{LinearSampling}, dist, survival) = cquantile(dist, survival)
invert_space(::Type{LogSampling}, dist, survival) = invlogccdf(dist, survival)


"""
This combines Next Reaction Method and Modified Next Reaction Method.
Each distribution is more precise being sampled in either a linear space
or a log space. This sampler chooses which space to use depending on the
type of the `UnivariateDistribution`. Defaults are set for those distributions
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
"""
struct CombinedNextReaction{T}
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    transition_entry::Dict{T,NRTransition}
end


function CombinedNextReaction{T}() where {T}
    heap = MutableBinaryMinHeap{OrderedSample{T}}()
    CombinedNextReaction{T}(heap, Dict{T,NRTransition}())
end


function next(nr::CombinedNextReaction{T}, when::Float64, rng::AbstractRNG) where {T}
    if !isempty(nr.firing_queue)
        least = top(nr.firing_queue)
        # For this sampler, mark this transition as the one that will fire
        # by marking its remaining cumulative time as 0.0.
        entry = nr.transition_entry[least.key]
        nr.transition_entry[least.key] = NRTransition(
            entry.heap_handle, get_survival_zero(nr), entry.distribution, entry.te, entry.t0
        )
        return (least.time, least.key)
    else
        # Return type is Tuple{Float64, Union{Nothing,T}} because T is not default-constructible.
        return (Inf, nothing)
    end
end


function sample_shifted(
    nr::CombinedNextReaction{T},
    rng::AbstractRNG,
    distribution::UnivariateDistribution,
    te::Float64,
    when::Float64
    ) where {T}
    if te < when
        shifted_distribution = truncated(distribution, when - te, Inf)
        sample = rand(rng, shifted_distribution)
        tau = te + sample
        survival = survival_space(shifted_distribution, sample)
    else  # te >= when
        # The distribution starts in the future
        sample = rand(rng, distribution)
        tau = te + sample
        survival = survival_space(distribution, sample)
    end
    (tau, survival)
end


function sample_by_inversion(
    nr::NextReaction{T},
    distribution::CombinedNextReaction, te::Float64, when::Float64, survival::Float64
    ) where {T}
    if te < when
        te + invert_space(truncated(distribution, when - te, Inf), survival)
    else   # te > when
        te + invert_space(distribution, survival)
    end
end


# Transition was enabled between time record.t0 and tn.
# Divide the survival by the conditional survival between t0 and tn.
# te can be before t0, at t0, between t0 and tn, or at tn, or after tn.
function consume_survival(nr::CombinedNextReaction{T}, record::NRTransition, tn::Float64) where {T}
    survive_te_tn = if record.te < tn
        survival_space(record.distribution, tn-record.te)
    else
        get_survival_max(record.distribution)
    end
    survive_te_t0 = if record.te < record.t0
        survival_space(record.distribution, record.t0-record.te)
    else
        get_survival_max(record.distribution)
    end
    record.survival / (survive_te_t0 * survive_te_tn)
end


function enable!(
    nr::CombinedNextReaction{T}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    # Three cases: a) never been enabled b) currently enabled c) was disabled.
    record = get(nr.transition_entry, clock, NRTransition(0, get_survival_zero(nr), Never(), 0.0, 0.0))
    heap_handle = record.heap_handle

    # if record.survival <= 0.0
    if record.survival <= get_survival_zero(nr)
        tau, survival = sample_shifted(nr, rng, distribution, te, when)
        sample = OrderedSample{T}(clock, tau)        
        if record.heap_handle > 0
            update!(nr.firing_queue, record.heap_handle, sample)
        else
            heap_handle = push!(nr.firing_queue, sample)
        end
        nr.transition_entry[clock] = NRTransition(
            heap_handle, survival, distribution, te, when
        )
    else
        # The transition was previously enabled.
        if record.heap_handle > 0
            # Consider te the same if the mantissa is within 2 bits of precision.
            same_te = abs(te - record.te) < 2 * eps(te)
            if same_te && distribution == record.distribution
                # No change. It's common to re-enable an already-enabled distribution.
            else
                # Account for time between when this was last enabled and now.
                survival = consume_survival(nr, record, when)
                tau = sample_by_inversion(nr, distribution, te, when, survival)
                entry = OrderedSample{T}(clock, tau)
                update!(nr.firing_queue, record.heap_handle, entry)
                nr.transition_entry[clock] = NRTransition(
                    heap_handle, survival, distribution, te, when
                )
            end

        # The transition was previously disabled.
        else
            tau = sample_by_inversion(nr, distribution, te, when, record.survival)
            heap_handle = push!(nr.firing_queue, OrderedSample{T}(clock, tau))
            nr.transition_entry[clock] = NRTransition(
                heap_handle, record.survival, distribution, te, when
            )
        end
    end
end


function disable!(nr::CombinedNextReaction{T}, clock::T, when::Float64) where {T}
    record = nr.transition_entry[clock]
    delete!(nr.firing_queue, record.heap_handle)
    nr.transition_entry[clock] = NRTransition(
        0, consume_survival(nr, record, when), record.distribution, record.te, when
    )
end
