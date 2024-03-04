import Base: sum!, push!, length, setindex!
using Random
using Distributions: Uniform
using Logging

"""
This is a binary tree where the leaves are values
and the nodes are sums of those values. It is meant to make it
easier to find the leaf such that the sum of it and all previous
leaves is greater than a given value.
"""
mutable struct BinaryTreePrefixSearch{T<:Real}
	array::Array{T,1}
	depth::Int64
	offset::Int64
	cnt::Int64
end


"""
    BinaryTreePrefixSearch{T}()

Constructor of a prefix search tree from an iterable list of real numbers.
"""
function BinaryTreePrefixSearch{T}(N=32) where {T<:Real}
	depth, offset, array_cnt = _btps_sizes(N)
	b = zeros(T, array_cnt)
	BinaryTreePrefixSearch{T}(b, depth, offset, 0)
end


time_type(ps::BinaryTreePrefixSearch{T}) where {T} = T
time_type(::Type{BinaryTreePrefixSearch{T}}) where {T} = T


function _btps_sizes(allocation)
	@assert allocation > 0
	depth = Int(ceil(log2(allocation))) + 1
	offset = 2^(depth - 1)
	array_cnt = 2^depth - 1
	@assert allocation <= offset + 1
	return (depth, offset, array_cnt)
end


# newsize is the desired number of entries. It will allocate more than this.
function resize!(pst::BinaryTreePrefixSearch{T}, newsize) where {T}
	depth, offset, array_cnt = _btps_sizes(newsize)
	b = zeros(T, array_cnt)
	will_fit = min(offset, pst.offset)
	b[offset:(offset + will_fit - 1)] = pst.array[pst.offset:(pst.offset + will_fit - 1)]
	pst.array = b
	pst.depth = depth
	pst.offset = offset
	pst.cnt = newsize
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
    for _ = 0:(pst.depth - 2)
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
	modify=Set{Int}()
	for (pos, value) in pairs
		index = pos + pst.offset - 1
		pst.array[index] = value
		push!(modify, div(index, 2))
	end

	# everything at depth-1 is correct, and changes are in modify.
    for depth = (pst.depth - 2):-1:0
        parents = Set{Int}()
        for node_idx in modify
            pst.array[node_idx] = pst.array[2 * node_idx] + pst.array[2 * node_idx + 1]
            push!(parents, div(node_idx, 2))
        end
        modify = parents
    end
end


function Base.push!(pst::BinaryTreePrefixSearch{T}, value::T) where T
	index = pst.cnt + 1
	if index <= allocated(pst)
		pst.cnt += 1
	else
		@debug "Pushing to binarytreeprefix $index $(allocated(pst))"
		resize!(pst, index)
	end
	pst[index] = value
	return value
end


"""
    setindex!(A, X, inds...)
"""
function Base.setindex!(pst::BinaryTreePrefixSearch{T}, value::T, index) where T
    set_multiple!(pst, [(index, value)])
end


# Private
function calculate_prefix!(pst::BinaryTreePrefixSearch)
	for d = (pst.depth - 2):-1:0
		for node_idx = (2^d):(2^(d + 1) - 1)
			pst.array[node_idx] = pst.array[2 * node_idx] + pst.array[2 * node_idx + 1]
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
