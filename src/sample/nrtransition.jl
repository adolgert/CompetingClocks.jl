# A record of a transition and the time.
# It's sortable by time. Immutable.
struct OrderedSample{K,T}
    key::K
    time::T
end

Base.:<(a::OrderedSample, b::OrderedSample) = a.time < b.time
Base.isless(a::OrderedSample, b::OrderedSample) = isless(a.time, b.time)
Base.:>(a::OrderedSample, b::OrderedSample) = a.time > b.time
Base.:(==)(a::OrderedSample, b::OrderedSample) = a.time == b.time
