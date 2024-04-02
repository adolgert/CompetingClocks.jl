using Random: rand, AbstractRNG
using Distributions: Uniform, Exponential, rate

export DirectCall, enable!, disable!, next


"""
    DirectCall{KeyType,TimeType,TreeType}()

DirectCall is responsible for sampling among Exponential distributions. It
samples using the Direct method. In this case, there is no optimization to
that Direct method, so we call it DirectCall because it recalculates
everything every time you call it.

The algorithm for the Direct Method relies heavily on what data structure
it uses to maintain a list of hazard rates, such that it can know the sum
of those hazards and index into them using a random value. This struct
has a default constructor that chooses a data structure for you, but there
are several options.

# Example

If we know that our simulation will only use a small number of different
clock keys, then it would make sense to use a data structure that disables
clocks by zeroing them out, instead of removing them from the list. This
will greatly reduce memory churn. We can do that by changing the
underlying data structure.

```julia
prefix_tree = BinaryTreePrefixSearch{T}()
keyed_prefix_tree = KeyedKeepPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
sampler_noremove = DirectCall{K,T,typeof(keyed_prefix_tree)}(keyed_prefix_tree)
```

"""
struct DirectCall{K,T,P} <: SSA{K,T}
    prefix_tree::P
end


function DirectCall{K,T}() where {K,T<:ContinuousTime}
    prefix_tree = BinaryTreePrefixSearch{T}()
    keyed_prefix_tree = KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    DirectCall{K,T,typeof(keyed_prefix_tree)}(keyed_prefix_tree)
end


reset!(dc::DirectCall) = (empty!(dc.prefix_tree); nothing)


"""
    enable!(dc::DirectCall, clock::T, distribution::Exponential, when, rng)

Tell the `DirectCall` sampler to enable this clock. The `clock` argument is
an identifier for the clock. The distribution is a univariate distribution
in time. In Julia, distributions are always relative to time `t=0`, but ours
start at some absolute enabling time, ``t_e``, so we provide that here.
The `when` argument is the time at which this clock is enabled, which may be
later than when it was first enabled. The `rng` is a random number generator.

If a particular clock had one rate before an event and it has another rate
after the event, call `enable!` to update the rate.
"""
function enable!(dc::DirectCall{K,T,P}, clock::K, distribution::Exponential,
        te::T, when::T, rng::AbstractRNG) where {K,T,P}
    dc.prefix_tree[clock] = rate(distribution)
end

function enable!(dc::DirectCall{K,T,P}, clock::K, distribution::D,
    te::T, when::T, rng::AbstractRNG) where {K,T,P,D<:UnivariateDistribution}
    error("DirectCall can only be used with Exponential type distributions")
end


"""
    disable!(dc::DirectCall, clock::T, when)

Tell the `DirectCall` sampler to disable this clock. The `clock` argument is
an identifier for the clock. The `when` argument is the time at which this
clock is enabled.
"""
function disable!(dc::DirectCall{K,T,P}, clock::K, when::T) where {K,T,P}
    delete!(dc.prefix_tree, clock)
end


"""
    next(dc::DirectCall, when::TimeType, rng::AbstractRNG)

Ask the sampler what clock will be the next to fire and at what time. This does
not change the sampler. You can call this multiple times and get multiple
answers. Each answer is a tuple of `(when, which clock)`. If there is no clock
to fire, then the response will be `(Inf, nothing)`. That's a good sign the
simulation is done.
"""
function next(dc::DirectCall{K,T,P}, when::T, rng::AbstractRNG) where {K,T,P}
    total = sum!(dc.prefix_tree)
    if total > eps(when)
        chosen, hazard_value = rand(rng, dc.prefix_tree)
        tau = when + rand(rng, Exponential(inv(total)))
        return (tau, chosen)
    else
        return (typemax(when), nothing)
    end
end

"""
    getindex(sampler::DirectCall{K,T}, clock::K)

For the `DirectCall` sampler, returns the rate parameter associated to the clock.
"""
function Base.getindex(dc::DirectCall{K,T,P}, clock::K) where {K,T,P}
    return dc.prefix_tree[clock]
end

function Base.keys(dc::DirectCall)
    return collect(keys(dc.prefix_tree.index))
end

function Base.length(dc::DirectCall)
    return length(dc.prefix_tree)
end
