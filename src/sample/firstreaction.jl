import Base: ==, <, >


"""
A record of a transition and the time.
It's sortable by time. Immutable.
"""
struct NRTransition{T}
	key::T
	time::Float64
end


function Base.:<(a::NRTransition, b::NRTransition)
	a.time < b.time
end


function Base.:>(a::NRTransition, b::NRTransition)
    a.time > b.time
end


function ==(a::NRTransition, b::NRTransition)
    a.time == b.time
end


"""
First reaction method for exponential transitions. At every time step, it
samples the distribution of every competing clock in order to find the soonest
to fire. This doesn't do any caching of rates. This sampler can be the fastest
when there are fewer than a dozen competing clocks.
"""
struct MarkovFirst{T}
    distributions::Dict{T,Float64}
end


function enable!(md::MarkovFirst{T}, clock::T, distribution::Exponential,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    md.distributions[clock] = params(distribution)[1]
end


function disable!(md::MarkovFirst{T}, clock::T, when::Float64) where {T}
    if haskey(md.distributions, clock)
        delete!(md.distributions, clock)
    end  # Else it didn't have this distribution listed.
end


"""
    next(md::MarkovFirst, process, when, rng::AbstractRNG)

This function is responsible for updating the state of the sampler and
returning both the next clock to fire and when it fires. The `process` is
the body of the simulation, its state and rules for what happens next.
The return value is `(time::Float64, clock_identifier)`.
That `process` must have a method called `hazards` which will report
all changes from the last event.

    hazards(process, rng::AbstractRNG, clock_update_function)

The `clock_update_function` takes the arguments
(`clock_identifier`, `distribution`, `enabled::Bool`).
"""
function next(md::MarkovFirst{T}, when::Float64, rng::AbstractRNG) where {T}
	soonest_clock = nothing
	soonest_time = Inf

    for (clock, propensity) in md.distributions
		# Don't use Gillespie's -log(rand()) sampling because it turns out
		# that's not what's used in practice.
		putative_time = rand(Exponential(propensity), rng)
		if putative_time < soonest_time
			soonest_clock = clock
			soonest_time = putative_time
		end
	end
	return (soonest_time, soonest_clock)
end


"""
Classic First Reaction method for Exponential and non-Exponential
distributions. Every time you sample, go to each distribution and ask when it
would fire. Then take the soonest and throw out the rest until the next sample.
"""
struct FirstReaction{T}
    distributions::Dict{T,(UnivariateDistribution,Float64)}
end


function enable!(md::FirstReaction{T}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {T}

    md.distributions[clock] = (distribution, te)
end


function disable!(md::FirstReaction{T}, clock::T, when::Float64) where {T}
    if haskey(md.distributions, clock)
        delete!(md.distributions, clock)
    end  # Else it didn't have this distribution listed.
end


function next(fr::FirstReaction, when::Float64, rng)
	soonest_clock = nothing
	soonest_time = Inf

    for (clock, (distribution, te)) in md.distributions
		relative_dist = truncated(distribution, when - te, Inf)
		putative_time = rand(relative_dist, rng)
		if putative_time < soonest_time
			soonest_clock = clock
			soonest_time = putative_time
		end
	end
	return (soonest_time, soonest_clock)
end
