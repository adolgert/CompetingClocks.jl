export SetOfSets

"""
A `SetOfSets` presents sets as a read-only unified set. The goal is to reduce memory allocation
by not presuming the client wants an instance of a Set object or an instance of
a Vector object. They may want to test is something is `in` the sets or want
to merge with another set.
"""
struct SetOfSets{K,T<:AbstractSet{K}} <: AbstractSet{K}
    subset::Vector{T}
    SetOfSets(subsets) = new{eltype(eltype(subsets)),eltype(subsets)}(collect(subsets))
end

function Base.iterate(sos::SetOfSets)
    isempty(sos.subset) && return nothing
    sub_idx = 1
    res = iterate(sos.subset[sub_idx])
    while res === nothing && sub_idx < lastindex(sos.subset)
        sub_idx += 1
        res = iterate(sos.subset[sub_idx])
    end
    res === nothing && return nothing
    return (res[1], (res[2], sub_idx))
end


function Base.iterate(sos::SetOfSets, (substate, sub_idx))
    res = iterate(sos.subset[sub_idx], substate)
    while res === nothing && sub_idx < lastindex(sos.subset)
        sub_idx += 1
        res = iterate(sos.subset[sub_idx])
    end
    res === nothing && return nothing
    return (res[1], (res[2], sub_idx))
end

Base.length(sos::SetOfSets) = sum(x -> length(x), sos.subset; init=0)
Base.in(x, sos::SetOfSets) = any(in(x, sub) for sub in sos.subset)
Base.eltype(::Type{SetOfSets{K,T}}) where {K,T} = K

Base.union(s::SetOfSets, other::SetOfSets) = SetOfSets(vcat(s.subset, other.subset))
Base.union(s::SetOfSets, other::AbstractSet) = SetOfSets(vcat(s.subset, [other]))
Base.intersect(s::SetOfSets, others...) = intersect(s.subset..., others...)
Base.setdiff(s::SetOfSets, others...) = setdiff(s.subset..., others...)
Base.issubset(s::SetOfSets, other) = all(issubset(x, other) for x in s.subset)

Base.isempty(sos::SetOfSets) = all(isempty, sos.subset)
