import Base: ==, <, >

"""
A record of a transition and the time.
It's sortable by time. Immutable.
"""
struct OrderedSample{T}
	key::T
	time::Float64
end


function Base.:<(a::OrderedSample, b::OrderedSample)
	a.time < b.time
end


function Base.isless(a::OrderedSample, b::OrderedSample)
    return isless(a.time, b.time)
end


function Base.:>(a::OrderedSample, b::OrderedSample)
    a.time > b.time
end


function ==(a::OrderedSample, b::OrderedSample)
    a.time == b.time
end
