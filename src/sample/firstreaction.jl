using Distributions: Exponential, params, truncated
using Random: rand, AbstractRNG

using Logging

export FirstReaction, ChatReaction

"""
    FirstReaction{KeyType,TimeType}()

This is the classic first reaction method for general distributions. 
Every time you sample, this goes to each distribution and asks
when it would fire. Then it takes the soonest and throws out the rest of the
sampled times until the next sample. It can also be very fast when there are only a few clocks to sample.

One interesting property of this sampler is that you can call `next()`
multiple times in order to get a distribution of next firing clocks and their
times to fire.
"""
struct FirstReaction{K,T} <: SSA{K,T}
	# This other class already stores the current set of distributions, so use it.
    core_matrix::TrackWatcher{K}
	FirstReaction{K,T}() where {K,T <: ContinuousTime} = new(TrackWatcher{K,T}())
end


reset!(fr::FirstReaction) = reset!(fr.core_matrix)
Base.copy!(dst::FirstReaction{K,T}, src::FirstReaction{K,T}) where {K,T} = (copy!(dst.core_matrix, src.core_matrix); dst)


function enable!(fr::FirstReaction{K,T}, clock::K, distribution::UnivariateDistribution,
		te::T, when::T, rng::AbstractRNG) where {K,T}

	enable!(fr.core_matrix, clock, distribution, te, when, rng)
end


function disable!(fr::FirstReaction{K,T}, clock::K, when::T) where {K,T}
	disable!(fr.core_matrix, clock, when)
end


function next(fr::FirstReaction{K,T}, when::T, rng::AbstractRNG) where {K,T}
	soonest_clock::Union{Nothing,K} = nothing
	soonest_time = typemax(T)

    for entry::EnablingEntry{K,T} in fr.core_matrix
		if entry.te < when
			relative_dist = truncated(entry.distribution, when - entry.te, typemax(T))
			putative_time = entry.te + rand(rng, relative_dist)
		else
			putative_time = entry.te + rand(rng, entry.distribution)
		end
		if putative_time < soonest_time
			soonest_clock = entry.clock
			soonest_time = putative_time
		end
	end
	return (soonest_time, soonest_clock)
end

"""
	getindex(sampler::FirstReaction{K,T}, clock::K)

For the `FirstReaction` sampler, returns the distribution object associated to the clock.
"""
function Base.getindex(fr::FirstReaction{K,T}, clock::K) where {K,T}
    if haskey(fr.core_matrix.enabled, clock)
        return getfield(fr.core_matrix.enabled[clock], :distribution)
    else
        throw(KeyError(clock))
    end
end

function Base.keys(fr::FirstReaction)
    return collect(keys(fr.core_matrix.enabled))
end

function Base.length(fr::FirstReaction)
    return length(fr.core_matrix.enabled)
end

Base.haskey(fr::FirstReaction, clock) = isenabled(fr.core_matrix, clock)


"""
This sampler can help if it's the first time you're trying a model. It checks
all of the things and uses Julia's logger to communicate them. It samples
using the first reaction algorithm.
"""
mutable struct ChatReaction{K,T}
	# This other class already stores the current set of distributions, so use it.
    core_matrix::TrackWatcher{K}
	step_cnt::Int64
	enables::Set{K}
	disables::Set{K}
	ChatReaction{K,T}() where {K,T<:ContinuousTime} = new(TrackWatcher{K,T}(), 0, Set{K}(), Set{K}())
end


function enable!(fr::ChatReaction{K,T}, clock::K, distribution::UnivariateDistribution,
		te::T, when::T, rng::AbstractRNG) where {K,T}

	if clock ∈ keys(fr.core_matrix.enabled)
		@warn "Re-enabling transition $clock without disabling first"
	end
	enable!(fr.core_matrix, clock, distribution, te, when, rng)
	push!(fr.enables, clock)
end


function disable!(fr::ChatReaction{K,T}, clock::K, when::T) where {K,T}
	if clock ∉ fr.enables
		@warn "Disabling a clock that was never enabled: $(clock)."
	end
	if clock ∉ keys(fr.core_matrix.enabled)
		@warn "Disabling clock $(clock) that wasn't currently enabled"
	end
	disable!(fr.core_matrix, clock, when)
	push!(fr.disables, clock)
end


function next(fr::ChatReaction{K,T}, when::T, rng) where {K,T}
	soonest_clock = nothing
	soonest_time = typemax(T)

	if length(fr.enables) == 0
		@warn "No transitions have ever been enabled. Sampler may not be initialized."
	end

    for entry::EnablingEntry{K,T} in fr.core_matrix
		if entry.te < when
			relative_dist = truncated(entry.distribution, when - entry.te, typemax(T))
			putative_time = entry.te + rand(rng, relative_dist)
		else
			putative_time = entry.te + rand(rng, entry.distribution)
		end
		if putative_time < soonest_time
			soonest_clock = entry.clock
			soonest_time = putative_time
		end
	end
	fr.step_cnt += 1
	@debug "Step $(fr.step_cnt) time $(soonest_time) fire $(soonest_clock)"
	return (soonest_time, soonest_clock)
end


Base.haskey(fr::ChatReaction, clock) = isenabled(fr.core_matrix, clock)
