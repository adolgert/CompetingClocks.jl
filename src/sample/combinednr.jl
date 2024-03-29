
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
function sampling_space(x::T) where {T <: UnivariateDistribution}
    sampling_space(typeof(x))::Union{Type{LinearSampling},Type{LogSampling}}
end
abstract type SamplingSpaceType end
struct LinearSampling <: SamplingSpaceType end
struct LogSampling <: SamplingSpaceType end
sampling_space(::Type) = LinearSampling
sampling_space(::Type{Distributions.Arcsine}) = LinearSampling
sampling_space(::Type{Distributions.BetaPrime}) = LinearSampling
sampling_space(::Type{Distributions.Biweight}) = LinearSampling
sampling_space(::Type{Distributions.Beta}) = LinearSampling
sampling_space(::Type{Distributions.Cauchy}) = LinearSampling
sampling_space(::Type{Distributions.Cosine}) = LinearSampling
sampling_space(::Type{Distributions.Epanechnikov}) = LinearSampling
sampling_space(::Type{Distributions.Erlang}) = LogSampling
sampling_space(::Type{Distributions.Exponential}) = LogSampling
sampling_space(::Type{Distributions.Frechet}) = LinearSampling
sampling_space(::Type{Distributions.Gamma}) = LogSampling
sampling_space(::Type{Distributions.GeneralizedPareto}) = LinearSampling
sampling_space(::Type{Distributions.Gumbel}) = LinearSampling
sampling_space(::Type{Distributions.InverseGamma}) = LinearSampling
sampling_space(::Type{Distributions.InverseGaussian}) = LinearSampling
sampling_space(::Type{Distributions.JohnsonSU}) = LinearSampling
sampling_space(::Type{Distributions.Kolmogorov}) = LinearSampling
sampling_space(::Type{Distributions.Kumaraswamy}) = LinearSampling
sampling_space(::Type{Distributions.Laplace}) = LogSampling
sampling_space(::Type{Distributions.Levy}) = LinearSampling
sampling_space(::Type{Distributions.Lindley}) = LinearSampling
sampling_space(::Type{Distributions.Logistic}) = LinearSampling
sampling_space(::Type{Distributions.LogitNormal}) = LinearSampling
sampling_space(::Type{Distributions.LogNormal}) = LinearSampling
sampling_space(::Type{Distributions.Normal}) = LinearSampling
sampling_space(::Type{Distributions.NormalCanon}) = LinearSampling
sampling_space(::Type{Distributions.Pareto}) = LinearSampling
sampling_space(::Type{Distributions.PGeneralizedGaussian}) = LinearSampling
sampling_space(::Type{Distributions.Rayleigh}) = LinearSampling
sampling_space(::Type{Distributions.Rician}) = LinearSampling
sampling_space(::Type{Distributions.SkewedExponentialPower}) = LinearSampling
sampling_space(::Type{Distributions.SymTriangularDist}) = LinearSampling
sampling_space(::Type{Distributions.Triweight}) = LinearSampling
sampling_space(::Type{Distributions.Uniform}) = LinearSampling
sampling_space(::Type{Distributions.Weibull}) = LogSampling

# The following four support functions are used by the CombinedNextReaction
# sampler, and their use is decided by the `sampling_space()` above.
# I had some trouble getting these functions to be type stable, and their lack
# of type stability hurt performance. Test performance by using
# tests/time_combinednr.jl.
function get_survival_zero(::T) where {T <: UnivariateDistribution}
    get_survival_zero(sampling_space(T))::Float64
end
function get_survival_zero(::Type{T}) where {T <: UnivariateDistribution}
    get_survival_zero(sampling_space(T))::Float64
end
get_survival_zero(::Type{LinearSampling}) = 0.0::Float64
get_survival_zero(::Type{LogSampling}) = -Inf::Float64

draw_space(::Type{LinearSampling}, rng) = rand(rng, Uniform())
draw_space(::Type{LogSampling}, rng) = rand(rng, Exponential())

function survival_space(::Type{T}, dist, sample) where {T <: UnivariateDistribution}
    survival_space(sampling_space(T), dist, sample)
end
survival_space(::Type{LinearSampling}, dist, sample) = ccdf(dist, sample)::Float64
survival_space(::Type{LogSampling}, dist, sample) = logccdf(dist, sample)::Float64

function invert_space(::Type{T}, dist, survival) where {T <: UnivariateDistribution}
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
    CombinedNextReaction{KeyType,TimeType}()

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
"""
struct CombinedNextReaction{K,T} <: SSA{K,T}
    firing_queue::MutableBinaryMinHeap{OrderedSample{K,T}}
    transition_entry::Dict{K,NRTransition{T}}
end


function CombinedNextReaction{K,T}() where {K,T <: ContinuousTime}
    heap = MutableBinaryMinHeap{OrderedSample{K,T}}()
    CombinedNextReaction{K,T}(heap, Dict{K,NRTransition{T}}())
end

clone(nr::CombinedNextReaction{K,T}) where {K,T} = CombinedNextReaction{K,T}()
export clone


@doc raw"""
For the first reaction sampler, you can call next() multiple times and get
different, valid, answers. That isn't the case here. When you call next()
on a CombinedNextReaction sampler, it returns the key associated with the
clock that fires and marks that clock as fired. Calling next() again would
return a nonsensical value.
"""
function next(nr::CombinedNextReaction{K,T}, when::T, rng::AbstractRNG) where {K,T}
    if !isempty(nr.firing_queue)
        least = first(nr.firing_queue)
        # For this sampler, mark this transition as the one that will fire
        # by marking its remaining cumulative time as 0.0.
        entry = nr.transition_entry[least.key]
        nr.transition_entry[least.key] = NRTransition{T}(
            entry.heap_handle, get_survival_zero(entry.distribution),
            entry.distribution, entry.te, entry.t0
        )
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
    ) where {S <: SamplingSpaceType, T <: ContinuousTime}
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
    ) where {S <: SamplingSpaceType, T <: ContinuousTime}
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
    ) where {S <: LinearSampling, T <: ContinuousTime}
    survive_te_tn = if record.te < tn
        ccdf(distribution, tn-record.te)::T
    else
        one(T)
    end
    survive_te_t0 = if record.te < record.t0
        ccdf(distribution, record.t0-record.te)::T
    else
        one(T)
    end
    record.survival / (survive_te_t0 * survive_te_tn)
end


@doc raw"""
This updates the survival for a transition in log space, according to
Anderson's method.

``\ln u=-\int_{t_e}^{t_n}\lambda_0(s-t_e)ds - \int_{t_n}^{\tau'}\lambda_{n}(s-t_e)ds``

"""
function consume_survival(
    record::NRTransition, distribution::UnivariateDistribution, ::Type{S}, tn::T
    ) where {S <: LogSampling, T <: ContinuousTime}
    log_survive_te_tn = if record.te < tn
        logccdf(distribution, tn-record.te)::T
    else
        zero(T)
    end
    log_survive_te_t0 = if record.te < record.t0
        logccdf(distribution, record.t0-record.te)::T
    else
        zero(T)
    end
    record.survival - (log_survive_te_t0 + log_survive_te_tn)
end


function enable!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T, rng::AbstractRNG) where {K,T}
    enable!(nr, clock, distribution, sampling_space(distribution), te, when, rng)
    nothing
end


function enable!(
    nr::CombinedNextReaction{K,T}, clock::K, distribution::UnivariateDistribution, ::Type{S},
    te::T, when::T, rng::AbstractRNG) where {K, T, S <: SamplingSpaceType}
    
    # Three cases: a) never been enabled b) currently enabled c) was disabled.
    record = get(
        nr.transition_entry,
        clock,
        NRTransition{T}(0, get_survival_zero(S), Never(), zero(T), zero(T))
        )
    heap_handle = record.heap_handle

    # if the transition needs to be re-drawn.
    if record.survival <= get_survival_zero(S)
        tau, shift_survival = sample_shifted(rng, distribution, S, te, when)
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


function disable!(nr::CombinedNextReaction{K,T}, clock::K, when::T) where {K,T <: ContinuousTime}
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
