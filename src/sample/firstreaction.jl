using Distributions: Exponential, params, truncated
using Random: rand, AbstractRNG

using Logging

export FirstReaction, ChatReaction

"""
Classic First Reaction method for Exponential and non-Exponential
distributions. Every time you sample, go to each distribution and ask when it
would fire. Then take the soonest and throw out the rest until the next sample.
"""
struct FirstReaction{T}
	# This other class already stores the current set of distributions, so use it.
    core_matrix::TrackWatcher{T}
	FirstReaction{T}() where {T} = new(TrackWatcher{T}())
end


function enable!(fr::FirstReaction{T}, clock::T, distribution::UnivariateDistribution,
		te::Float64, when::Float64, rng::AbstractRNG) where {T}

	enable!(fr.core_matrix, clock, distribution, te, when, rng)
end


function disable!(fr::FirstReaction{T}, clock::T, when::Float64) where {T}
	disable!(fr.core_matrix, clock, when)
end


function next(fr::FirstReaction{T}, when::Float64, rng) where {T}
	soonest_clock = nothing
	soonest_time = Inf

    for entry::EnablingEntry{T} in fr.core_matrix
		if entry.te < when
			relative_dist = truncated(entry.distribution, when - entry.te, Inf)
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
This sampler can help if it's the first time you're trying a model. It checks
all of the things and uses Julia's logger to communicate them.
"""
mutable struct ChatReaction{T}
	# This other class already stores the current set of distributions, so use it.
    core_matrix::TrackWatcher{T}
	step_cnt::Int64
	enables::Set{T}
	disables::Set{T}
	ChatReaction{T}() where {T} = new(TrackWatcher{T}(), 0, Set{T}(), Set{T}())
end


function enable!(fr::ChatReaction{T}, clock::T, distribution::UnivariateDistribution,
		te::Float64, when::Float64, rng::AbstractRNG) where {T}

	if clock ∈ keys(fr.core_matrix.enabled)
		@warn "Re-enabling transition $clock without disabling first"
	end
	enable!(fr.core_matrix, clock, distribution, te, when, rng)
	push!(fr.enables, clock)
end


function disable!(fr::ChatReaction{T}, clock::T, when::Float64) where {T}
	if clock ∉ fr.enables
		@warn "Disabling a clock that was never enabled: $(clock)."
	end
	if clock ∉ keys(fr.core_matrix.enabled)
		@warn "Disabling clock $(clock) that wasn't currently enabled"
	end
	disable!(fr.core_matrix, clock, when)
	push!(fr.disables, clock)
end


function next(fr::ChatReaction{T}, when::Float64, rng) where {T}
	soonest_clock = nothing
	soonest_time = Inf

	if length(fr.enables) == 0
		@warn "No transitions have ever been enabled. Sampler may not be initialized."
	end

    for entry::EnablingEntry{T} in fr.core_matrix
		if entry.te < when
			relative_dist = truncated(entry.distribution, when - entry.te, Inf)
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
