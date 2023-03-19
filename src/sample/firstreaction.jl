using Distributions: Exponential, params, truncated
using Random: rand, AbstractRNG

export FirstReaction

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


function disable!(md::FirstReaction{T}, clock::T, when::Float64) where {T}
	disable!(md.core_matrix, clock, when)
end


function next(fr::FirstReaction{T}, when::Float64, rng) where {T}
	soonest_clock = nothing
	soonest_time = Inf

    for entry::EnablingEntry{T} in fr.core_matrix
		relative_dist = truncated(entry.distribution, when - entry.te, Inf)
		putative_time = rand(rng, relative_dist)
		if putative_time < soonest_time
			soonest_clock = entry.clock
			soonest_time = putative_time
		end
	end
	return (soonest_time, soonest_clock)
end
