# A record of a transition and the time.
# It's sortable by time. Immutable.
struct OrderedSample{K,T}
    key::K
    time::T
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


function Base.:(==)(a::OrderedSample, b::OrderedSample)
    a.time == b.time
end
