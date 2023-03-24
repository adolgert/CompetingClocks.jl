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


function Base.isless(a::NRTransition, b::NRTransition)
    return isless(a.time, b.time)
end


function Base.:>(a::NRTransition, b::NRTransition)
    a.time > b.time
end


function ==(a::NRTransition, b::NRTransition)
    a.time == b.time
end
