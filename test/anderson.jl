using DataStructures
using Random

using Fleck: OrderedSample


"""
    MNRM by Anderson
"""
struct ModifiedNextReaction{T} <: AbstractNextReaction{T}
    firing_queue::MutableBinaryHeap{OrderedSample{T}}
    transition_entry::Dict{T,NRTransition}
end

get_survival_zero(::ModifiedNextReaction{T}) where {T} = -Inf

function ModifiedNextReaction{T}() where {T}
    heap = MutableBinaryMinHeap{OrderedSample{T}}()
    ModifiedNextReaction{T}(heap, Dict{T,NRTransition}())
end

# next (use abstract method)


function sample_shifted(
    nr::ModifiedNextReaction{T},
    rng::AbstractRNG, distribution::UnivariateDistribution, te::Float64, when::Float64
    ) where {T}
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

function sample_by_inversion(
    nr::ModifiedNextReaction{T},
    distribution::UnivariateDistribution, te::Float64, when::Float64, logsurvival::Float64
    ) where {T}
    if te < when
        te + invlogccdf(truncated(distribution, when - te, Inf), logsurvival)
    else   # te > when
        te + invlogccdf(distribution, logsurvival)
    end
end

function consume_survival(nr::ModifiedNextReaction{T}, record::NRTransition, tn::Float64) where {T}
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
    record.survival - (log_survive_te_t0 + log_survive_te_tn)
end

# const MNRNotFound = MNRTransition(0, -Inf, Never(), 0.0, 0.0)

# function enable!(
#     nr::ModifiedNextReaction{T}, clock::T, distribution::UnivariateDistribution,
#     te::Float64, when::Float64, rng::AbstractRNG) where {T}

#     # Three cases: a) never been enabled b) currently enabled c) was disabled.
#     record = get(nr.transition_entry, clock, MNRNotFound)
#     heap_handle = record.heap_handle

#     if record.log_survival <= -Inf
#         tau, log_survival = sample_shifted(nr, rng, distribution, te, when)
#         sample = OrderedSample{T}(clock, tau)        
#         if record.heap_handle > 0
#             update!(nr.firing_queue, record.heap_handle, sample)
#         else
#             heap_handle = push!(nr.firing_queue, sample)
#         end
#         nr.transition_entry[clock] = MNRTransition(
#             heap_handle, log_survival, distribution, te, when
#         )
#     else
#         # The transition was previously enabled.
#         if record.heap_handle > 0
#             # Consider te the same if the mantissa is within 2 bits of precision.
#             same_te = abs(te - record.te) < 2 * eps(te)
#             if same_te && distribution == record.distribution
#                 # No change. It's common to re-enable an already-enabled distribution.
#             else
#                 # Account for time between when this was last enabled and now.
#                 log_survival = consume_survival(record, when)
#                 tau = sample_by_inversion(nr, distribution, te, when, log_survival)
#                 entry = OrderedSample{T}(clock, tau)
#                 update!(nr.firing_queue, record.heap_handle, entry)
#                 nr.transition_entry[clock] = MNRTransition(
#                     heap_handle, log_survival, distribution, te, when
#                 )
#             end
    

#         # The transition was previously disabled.
#         else
#             tau = sample_by_inversion(nr, distribution, te, when, record.log_survival)
#             heap_handle = push!(nr.firing_queue, OrderedSample{T}(clock, tau))
#             nr.transition_entry[clock] = MNRTransition(
#                 heap_handle, record.log_survival, distribution, te, when
#             )
#         end
#     end
# end

# function disable!(nr::ModifiedNextReaction{T}, clock::T, when::Float64) where {T}
#     record = nr.transition_entry[clock]
#     delete!(nr.firing_queue, record.heap_handle)
#     nr.transition_entry[clock] = MNRTransition(
#         0, consume_survival(record, when), record.distribution, record.te, when
#     )
# end