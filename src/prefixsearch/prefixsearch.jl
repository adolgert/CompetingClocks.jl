import Base: sum, push!

"""
This is a binary tree where the leaves are values
and the nodes are sums of those values. It is meant to make it
easier to find the leaf such that the sum of it and all previous
leaves is greater than a given value.
"""
struct PrefixSearchTree{T<:Real}
	array::Array{T,1}
	depth::Int
	offset::Int
end


"""
    PrefixSearchTree(ElementType::DataType, count::Int)

Constructor for an empty tree of type ElementType and length count.
"""
function PrefixSearchTree(ElementType::DataType, count::Int)
	new(zeros(ElementType, count), 0, 0)
end


"""
    PrefixSearchTree(summables)

Constructor of a prefix search tree from an iterable list of real numbers.
"""
function PrefixSearchTree(summables)
	a = Array(summables)
	@assert length(a) > 0
	depth = Int(ceil(log(2, length(a)))) + 1
	offset = 2^(depth - 1)
	b = zeros(eltype(summables), 2^depth - 1)
	b[offset:(length(a) + offset - 1)] = a
	#print("depth $depth offset $offset length $(length(b))\n")
	pst = PrefixSearchTree{eltype(summables)}(b, depth, offset)
	calculate_prefix!(pst)
	pst
end


"""
    choose(pst::PrefixSearchTree, value)

Find the minimum index such that the prefix is greater than the given
value.

Precondition: The value must be strictly less than the total for
the tree.
"""
function choose(pst::PrefixSearchTree, value)
	if value >= pst.array[1]
		error("Value $value not less than total $(pst.array[1])")
	end

	requested_value = value
    index = 1
    for level = 0:(pst.depth - 2)
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


sum(pst::PrefixSearchTree) = pst.array[1]


"""
If there are multiple values to enter, then present them
at once as pairs of tuples, (index, value).
"""
function push!(pst::PrefixSearchTree, pairs)
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
            pst.array[node_idx]=pst.array[2 * node_idx] + pst.array[2 * node_idx + 1]
            push!(parents, div(node_idx, 2))
        end
        modify = parents
    end
end


function push!(pst::PrefixSearchTree{T}, index::Int, value::T) where T
    push!(pst, [(index, value)])
end


# Private
function calculate_prefix!(pst::PrefixSearchTree)
	for d = (pst.depth - 2):-1:0
		for node_idx = (2^d):(2^(d + 1) - 1)
			pst.array[node_idx] = pst.array[2 * node_idx] + pst.array[2 * node_idx + 1]
		end
	end
end
