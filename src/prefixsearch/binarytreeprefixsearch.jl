import Base: sum!, push!, length, setindex!, getindex
using Random
using Distributions: Uniform
using Logging

"""
    BinaryTreePrefixSearch{T}(N=32)

This stores hazard rates to make them faster for the Direct method to sample.
This is a binary tree where the leaves are values
and the nodes are sums of those values. It is meant to make it
easier to find the leaf such that the sum of it and all previous
leaves is greater than a given value.
"""
mutable struct BinaryTreePrefixSearch{T<:Real}
    # Data structure uses an array with children of i at 2i and 2i+1.
    array::Array{T,1} # length(array) > 0
    depth::Int # length(array) = 2^depth - 1
    offset::Int # 2^(depth - 1). Index of first leaf and number of leaves.
    cnt::Int # Number of leaves in use. Logical number of entries. cnt > 0.
    initial_allocation::Int
end


"""
    BinaryTreePrefixSearch{T}([N])

Constructor of a prefix search tree from an iterable list of real numbers.
The optional hint, N, is the number of values to pre-allocate.
"""
function BinaryTreePrefixSearch{T}(N=32) where {T<:Real}
    depth, offset, array_cnt = _btps_sizes(N)
    b = zeros(T, array_cnt)
    BinaryTreePrefixSearch{T}(b, depth, offset, 0, N)
end


function Base.empty!(ps::BinaryTreePrefixSearch)
    depth, offset, array_cnt = _btps_sizes(ps.initial_allocation)
    resize!(ps.array, array_cnt)
    fill!(ps.array, zero(eltype(ps.array)))
    ps.depth = depth
    ps.offset = offset
    ps.cnt = 0
end

function Base.copy!(dst::BinaryTreePrefixSearch{T}, src::BinaryTreePrefixSearch{T}) where {T}
    copy!(dst.array, src.array)
    dst.depth = src.depth
    dst.offset = src.offset
    dst.cnt = src.cnt
    dst.initial_allocation = src.initial_allocation
    return dst
end


time_type(ps::BinaryTreePrefixSearch{T}) where {T} = T
time_type(::Type{BinaryTreePrefixSearch{T}}) where {T} = T


# Integer log2, rounding up.
function ceil_log2(v::Integer)
    r = 0
    power_of_two = ((v & (v - 1)) == 0) ? 0 : 1
    while (v >>= 1) != 0
        r += 1
    end
    r + power_of_two
end


# The tree must have at least `allocation` leaves.
function _btps_sizes(allocation)
    allocation = (allocation > 0) ? allocation : 1
    depth = ceil_log2(allocation) + 1
    offset = 2^(depth - 1)
    array_cnt = 2^depth - 1
    @assert allocation <= offset + 1
    return (depth, offset, array_cnt)
end


# newcnt is the desired number of entries.
function Base.resize!(pst::BinaryTreePrefixSearch{T}, newcnt) where {T}
    depth, offset, array_cnt = _btps_sizes(newcnt)
    b = zeros(T, array_cnt)
    will_fit = min(offset, pst.offset)
    b[offset:(offset+will_fit-1)] = pst.array[pst.offset:(pst.offset+will_fit-1)]
    pst.array = b
    pst.depth = depth
    pst.offset = offset
    pst.cnt = newcnt
    calculate_prefix!(pst)
end


Base.length(ps::BinaryTreePrefixSearch) = ps.cnt
allocated(ps::BinaryTreePrefixSearch) = ps.offset


"""
    choose(pst::BinaryTreePrefixSearch, value)

Find the minimum index such that the prefix is greater than the given
value.

Precondition: The value must be strictly less than the total for
the tree.
"""
function choose(pst::BinaryTreePrefixSearch{T}, value) where {T}
    if value >= pst.array[1]
        error("Value $value not less than total $(pst.array[1])")
    end

    index = 1
    for _ = 0:(pst.depth-2)
        left_child = 2 * index
        if pst.array[left_child] > value
            index = left_child
        else
            index = left_child + 1
            value -= pst.array[left_child]
        end
    end
    index - pst.offset + 1, pst.array[index] - value
end


Base.sum!(pst::BinaryTreePrefixSearch) = pst.array[1]


"""
If there are multiple values to enter, then present them
at once as pairs of tuples, (index, value).
"""
function set_multiple!(pst::BinaryTreePrefixSearch, pairs)
    maxindex = maximum([i for (i, v) in pairs])
    if maxindex > allocated(pst)
        @debug "BinaryTreePrefixSearch resizing to $maxindex"
        resize!(pst, maxindex)
    end
    if maxindex > pst.cnt
        pst.cnt = maxindex
    end
    modify = Set{Int}()
    sizehint!(modify, length(pairs))
    for (pos, value) in pairs
        index = pos + pst.offset - 1
        pst.array[index] = value
        push!(modify, div(index, 2))
    end

    for depth = (pst.depth-1):-1:1
        parents = Set{Int}()
        for node_idx in modify
            pst.array[node_idx] = pst.array[2*node_idx] + pst.array[2*node_idx+1]
            push!(parents, div(node_idx, 2))
        end
        modify = parents
    end
    @assert length(modify) == 1 && first(modify) == 0
end


function Base.push!(pst::BinaryTreePrefixSearch{T}, value::T) where T
    set_multiple!(pst, [(pst.cnt + 1, value)])
    return value
end


"""
    setindex!(A, X, inds...)
"""
function Base.setindex!(pst::BinaryTreePrefixSearch{T}, value::T, index) where T
    set_multiple!(pst, [(index, value)])
end

function Base.getindex(pst::BinaryTreePrefixSearch{T}, index) where {T}
    return pst.array[index+pst.offset-1]
end


# Private
function calculate_prefix!(pst::BinaryTreePrefixSearch)
    for d = (pst.depth-2):-1:0
        for node_idx = (2^d):(2^(d+1)-1)
            pst.array[node_idx] = pst.array[2*node_idx] + pst.array[2*node_idx+1]
        end
    end
end


"""
    rand(rng, sampler::SamplerTrivial{BinaryTreePrefixSearch})

This method overload allows the machinery of Random to generate random variates
from the BinaryTreePrefixSearch set of values.
"""
Random.rand(rng::AbstractRNG, d::Random.SamplerTrivial{BinaryTreePrefixSearch{T}}) where {T} =
    choose(d[], rand(rng, Uniform{T}(zero(T), d[].array[1])))


Base.haskey(md::BinaryTreePrefixSearch, clock) = false

function Base.haskey(pst::BinaryTreePrefixSearch{T}, clock::Int) where {T}
    if 0 < clock ≤ length(pst)
        return getindex(pst, clock) > zero(T)
    else
        return false
    end
end
