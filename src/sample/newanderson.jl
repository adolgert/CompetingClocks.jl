using DataStructures

export ModifiedNextReaction

struct MNRTransition
    heap_handle::Int
    log_survival::Float64 # when cumulative intensity = log surv, it fires
    distribution::UnivariateDistribution
    te::Float64  # Enabling time of distribution
    t0::Float64  # Enabling time of transition
end


"""
    MNRM by Anderson
"""
struct ModifiedNextReaction{T}
    # This stores everything we need to know about enabled transitions.
    # This struct requires the handle of the heap to always be greater than zero.
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    # This maps from transition to entry in the firing queue and remaining
    # integrated hazard. Has info on _disabled_ transitions.
    transition_entry::Dict{T,MNRTransition}
end

function ModifiedNextReaction{T}() where {T}
    heap = MutableBinaryMinHeap{OrderedSample{T}}()
    ModifiedNextReaction{T}(heap, Dict{T,MNRTransition}())
end

function next(nr::ModifiedNextReaction{T}, when::Float64, rng::AbstractRNG) where {T}
    if !isempty(nr.firing_queue)
        least = top(nr.firing_queue)
        entry = nr.transition_entry[least.key]
        nr.transition_entry[least.key] = MNRTransition(
            # entry.heap_handle, 0.0, entry.distribution, entry.te, entry.t0
            entry.heap_handle, -Inf, entry.distribution, entry.te, entry.t0
        )
        return (least.time, least.key)
    else
        # Return type is Tuple{Float64, Union{Nothing,T}} because T is not default-constructible.
        return (Inf, nothing)
    end
end


function sample_shifted_mnrm(
    rng::AbstractRNG, distribution::UnivariateDistribution, te::Float64, when::Float64
    )
    if te < when
        shifted_distribution = truncated(distribution, when - te, Inf)
        sample = rand(rng, shifted_distribution)
        tau = te + sample
        log_survival = logccdf(shifted_distribution, sample)
    else  # te >= when
        # The distribution starts in the future
        sample = rand(rng, distribution)
        tau = te + sample
        log_survival = logccdf(distribution, sample)
    end
    (tau, log_survival)
end

function sample_by_inversion_mnrm(
    distribution::UnivariateDistribution, te::Float64, when::Float64, logsurvival::Float64
    )
    if te < when
        te + invlogccdf(truncated(distribution, when - te, Inf), logsurvival)
    else   # te > when
        te + invlogccdf(distribution, logsurvival)
    end
end

function consume_log_survival(record::MNRTransition, tn::Float64)
    log_survive_te_tn = if record.te < tn
        logccdf(record.distribution, tn-record.te)
    else
        0
    end
    log_survive_te_t0 = if record.te < record.t0
        logccdf(record.distribution, record.t0-record.te)
    else
        0
    end
    record.log_survival - (log_survive_te_t0 + log_survive_te_tn)
end

# const MNRNotFound = MNRTransition(0, -1.0, Never(), 0.0, 0.0)
const MNRNotFound = MNRTransition(0, -Inf, Never(), 0.0, 0.0)

function enable!(
    nr::ModifiedNextReaction{T}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    # Three cases: a) never been enabled b) currently enabled c) was disabled.
    record = get(nr.transition_entry, clock, MNRNotFound)
    heap_handle = record.heap_handle

    # if record.log_survival <= 0
    if record.log_survival <= -Inf
        tau, log_survival = sample_shifted_mnrm(rng, distribution, te, when)
        sample = OrderedSample{T}(clock, tau)        
        if record.heap_handle > 0
            update!(nr.firing_queue, record.heap_handle, sample)
        else
            heap_handle = push!(nr.firing_queue, sample)
        end
        nr.transition_entry[clock] = MNRTransition(
            heap_handle, log_survival, distribution, te, when
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
                log_survival = consume_log_survival(record, when)
                tau = sample_by_inversion_mnrm(distribution, te, when, log_survival)
                entry = OrderedSample{T}(clock, tau)
                update!(nr.firing_queue, record.heap_handle, entry)
                nr.transition_entry[clock] = MNRTransition(
                    heap_handle, log_survival, distribution, te, when
                )
            end

        # The transition was previously disabled.
        else
            tau = sample_by_inversion_mnrm(distribution, te, when, record.log_survival)
            heap_handle = push!(nr.firing_queue, OrderedSample{T}(clock, tau))
            nr.transition_entry[clock] = MNRTransition(
                heap_handle, record.log_survival, distribution, te, when
            )
        end
    end
end

function disable!(nr::ModifiedNextReaction{T}, clock::T, when::Float64) where {T}
    record = nr.transition_entry[clock]
    delete!(nr.firing_queue, record.heap_handle)
    nr.transition_entry[clock] = MNRTransition(
        0, consume_log_survival(record, when), record.distribution, record.te, when
    )
end



# test(::T) where {T} = @error("test not implemented for type $(T)")
# test(::Int) = println("you gave me an int!")

# test("hi")
# test(5)

# distribution = Weibull(2,5.0)

# # how to sample by inversion
# U = rand() # U
# samp = quantile(distribution,U) # F^{-1}(x) |-> y ∈ Ω

# log_s = -log(1-U) # 1-U = survival, -log(1-U) = samp log survival
# invlogccdf(distribution, -log_s)


# samp = rand(distribution)
# invlogccdf(distribution, logccdf(distribution, samp))
