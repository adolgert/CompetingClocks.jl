# Prefix Search

## What is it

Prefix search is also known as prefix sum, cumulative sum, or inclusive scan. These prefix search algorithms are responsible for

 * Calculating the sum of a list of values
 * Mutating values, where a mutated value may have a different index
   within the list
 * Finding the index of a value by its cumulant within the list

In Fleck, the Direct method uses prefix search to find the next hazard, or propensity, within the list of all enabled hazards. Most variations on the Direct method are optimizations of the prefix search algorithm.

## Example

For an example of use, let's create two propensities and select one.
```julia
ps = PrefixSearch()
push!(ps, propensity1)
push!(ps, propensity2)
uniform_variate1 = rand(rng, Uniform(0, sum(ps)))
chosen1 = choose(ps, uniform_variate1)
ps[1] = propensity3
uniform_variate2 = rand(rng, Uniform(0, sum(ps)))
chosen2 = choose(ps, uniform_variate2)
```

## Points of Variation

For a Direct method within a sampler framework, we are looking for a few features in our prefix search algorithm.

 * It may be indexed with a key rather than an integer.
 * There may be a fixed number of values, or those values may grow and shrink.
 * The size of values may be very small or very large, affecting roundoff error.

### Non-integer keys

I don't know of any implementations that can use non-integer keys, but Fleck allows any kind of clock ID you want, not just integers. A prefix search tree for non-integer keys would be some combination of a hash, or dict, and a set of cumulative sums. Julia does have sorted dictionaries in its DataStructures.jl library, but these are sorted by the key, not the value, so it's backwards from what we need.

### Mutable values

Some implementations never remove values. Instead, they set the propensity to zero and leave an entry in the table of propensities. Other implementations do remove the value or do it on a schedule, like rebalancing a tree data structure. Which method is fastest depends on the simulation, so it would be good to allow both.

Let's say there are two types of prefix search data structures, those that can remove an item and those that can only set it to zero. For those that can only set it to zero, the Direct method will keep track of what was removed and reuse those indices. For those that can remove items, they will support a `deleteat!()` method.


### Different sizes of values

If some propensities are very small and others very large, then adding the small to the large won't result in a change to the cumulant. This is a worse problem for 32-bit floating-point, so I don't think it's a large risk now for Float64, but historical solutions included sorting propensities, so that small is added to small (as in SortedDirect), or using tree structures that use multiple additions in order to mitigate the problem (as in OptimizedDirect).


## Constructing an Interface

I want to keep the core algorithms clear, so let's assemble this in layers.

 1. PrefixSearch - Whatever data structure underlies this, it looks to the client like a list of hazards. That means it has an integer index, from 1 to its length, and it can grow. We won't implement shrinking at this point.

    * setindex! - Set value of a hazard
    * push! - Add a new hazard
    * length - Total number of hazards
    * choose - Select a hazard by its cumulant
    * sum! - Find total of all hazards

 2. KeyedKeepPrefixSearch - To the client, this looks like a dictionary of hazards. It has an arbitrary type for the index. New keys can be added and keys can be deleted.

    * setindex! - Set value of a hazard
    * delete! - Disable a hazard
    * choose - Select a hazard by its cumulant
    * sum! - Find total of all hazards
