using DataStructures

export NextReaction

abstract type AbstractNextReaction{T} end


struct NRTransition
    heap_handle::Int
    survival::Float64 # value of S_{j}
    distribution::UnivariateDistribution
    te::Float64  # Enabling time of distribution
    t0::Float64  # Enabling time of transition
end

"""
    NextReaction{KeyType}()

This sampler uses the Next Reaction method popularized by Gibson and Bruck
in their 2000 paper called ``Efficient Exact Stochastic Simulation of Chemical
Systems with Many Species and Many Channels." The theory behind this method is
that drawing random numbers is the slowest part of sampling for the next time,
so they will minimize drawing random numbers. Currently random number generation
isn't the slowest part. It's inversion of the distributions, and there is a
lot of that with this method.
"""
struct NextReaction{T} <: AbstractNextReaction{T}
    # This stores everything we need to know about enabled transitions.
    # This struct requires the handle of the heap to always be greater than zero.
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    # This maps from transition to entry in the firing queue and remaining
    # integrated hazard. Has info on _disabled_ transitions.
    transition_entry::Dict{T,NRTransition}
end

get_survival_zero(::NextReaction{T}) where {T} = 0.0


function NextReaction{T}() where {T}
    heap = MutableBinaryMinHeap{OrderedSample{T}}()
    NextReaction{T}(heap, Dict{T,NRTransition}())
end


function next(nr::AbstractNextReaction{T}, when::Float64, rng::AbstractRNG) where {T}
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


"""
    sample_shifted(rng, distribution, te::Float64, when::Float64)

Sample from a distribution where the distribution and the sampling happen at
some absolute time. The distribution is specified as a distribution that starts
at a traditional zero value and a `te` that is the absolute time of that zero.
The `when` time is when we sample the distribution, going forward. Usually, `te`
is in the past, but it could be greater than `when`.

This samples a distribution conditional on the distribution not having fired
between time 0 and time ``(when - te)``, for ``te < when``

Returns the absolute time at which this distribution would fire.
"""
function sample_shifted(
    nr::NextReaction{T},
    rng::AbstractRNG, distribution::UnivariateDistribution, te::Float64, when::Float64
    ) where {T}
    if te < when
        shifted_distribution = truncated(distribution, when - te, Inf)
        sample = rand(rng, shifted_distribution)
        tau = te + sample
        survival = ccdf(shifted_distribution, sample)
    else  # te >= when
        # The distribution starts in the future
        sample = rand(rng, distribution)
        tau = te + sample
        survival = ccdf(distribution, sample)
    end
    (tau, survival)
end


function sample_by_inversion(
    nr::NextReaction{T},
    distribution::UnivariateDistribution, te::Float64, when::Float64, survival::Float64
    ) where {T}
    if te < when
        te + cquantile(truncated(distribution, when - te, Inf), survival)
    else   # te > when
        te + cquantile(distribution, survival)
    end
end


# Transition was enabled between time record.t0 and tn.
# Divide the survival by the conditional survival between t0 and tn.
# te can be before t0, at t0, between t0 and tn, or at tn, or after tn.
function consume_survival(nr::NextReaction{T}, record::NRTransition, tn::Float64) where {T}
    survive_te_tn = if record.te < tn
        ccdf(record.distribution, tn-record.te)
    else
        1
    end
    survive_te_t0 = if record.te < record.t0
        ccdf(record.distribution, record.t0-record.te)
    else
        1
    end
    record.survival / (survive_te_t0 * survive_te_tn)
end


# const NRNotFound = NRTransition(0, -1.0, Never(), 0.0, 0.0)


function enable!(
    nr::AbstractNextReaction{T}, clock::T, distribution::UnivariateDistribution,
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


function disable!(nr::AbstractNextReaction{T}, clock::T, when::Float64) where {T}
    record = nr.transition_entry[clock]
    delete!(nr.firing_queue, record.heap_handle)
    nr.transition_entry[clock] = NRTransition(
        0, consume_survival(nr, record, when), record.distribution, record.te, when
    )
end
