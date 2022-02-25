import Base: ==, <, >


"""
A record of a transition and the time.
It's sortable by time. Immutable.
"""
struct NRTransition
	key
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
Classic First Reaction method
"""
struct FirstReaction
end


function Next(fr::FirstReaction, system, rng)
	least=NRTransition(nothing, Inf)
	Hazards(system, rng) do clock, now, enabled, rng2
	  trial_time = Sample(clock.intensity, now, rng2)
	  @assert(trial_time >= now)
	  if trial_time < least.time
	  	least = NRTransition(clock, trial_time)
	  end
    end
    (least.time, least.key)
end


Observer(fr::FirstReaction) = (hazard, time, updated, rng) -> nothing
