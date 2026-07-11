using Random: rand, AbstractRNG
using Distributions: Uniform, Exponential, rate

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
mutable struct DirectCall{K,T,P} <: SSA{K,T}
    prefix_tree::P
    now::T
    log_likelihood::Float64
    calculate_likelihood::Bool
    streams::KeyedStreams{K}
end


function DirectCall{K,T}(; trajectory=false, seed=_DEFAULT_STREAM_SEED, coupling::Symbol=:redraw) where {K,T<:ContinuousTime}
    # Validated but not stored: a memoryless sampler can only redraw, so the
    # keyword exists to reject coupling=:carry at construction.
    validate_coupling(DirectCall, coupling)
    prefix_tree = BinaryTreePrefixSearch{T}()
    keyed_prefix_tree = KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    DirectCall{K,T,typeof(keyed_prefix_tree)}(keyed_prefix_tree, 0.0, 0.0, trajectory, KeyedStreams{K}(seed))
end


function DirectCallExplicit(
    ::Type{K}, ::Type{T}, ::Type{Keep}, ::Type{Prefix};
    trajectory=false, seed=_DEFAULT_STREAM_SEED, coupling::Symbol=:redraw) where {K,T,Keep,Prefix}

    validate_coupling(DirectCall, coupling)
    prefix_tree = Prefix{T}()
    keyed_prefix_tree = Keep{K,typeof(prefix_tree)}(prefix_tree)
    DirectCall{K,T,typeof(keyed_prefix_tree)}(keyed_prefix_tree, 0.0, 0.0, trajectory, KeyedStreams{K}(seed))
end


# DirectCall keeps no per-clock draw — it draws the next time and the winner from
# the race stream at each `next` — so a full-state clone is the clock table plus
# the race generator's state. The clone continues the same race draws.
function clone(dc::DirectCall{K,T,P}) where {K,T,P}
    c = DirectCall{K,T}(trajectory=dc.calculate_likelihood, seed=dc.streams.seed)
    copy!(c.prefix_tree, dc.prefix_tree)
    c.now = dc.now
    c.log_likelihood = dc.log_likelihood
    c.streams = copy(dc.streams)
    return c
end

similar_sampler(dc::DirectCall{K,T,P}) where {K,T,P} =
    DirectCall{K,T}(trajectory=dc.calculate_likelihood, seed=dc.streams.seed)

rekey_streams!(dc::DirectCall, seed) = (rekey_streams!(dc.streams, seed); dc)


# Nothing to resample: DirectCall draws afresh from the race stream at each next.
jitter!(dc::DirectCall{K,T,P}, when::T) where {K,T,P} = nothing


function reset!(dc::DirectCall{K,T,P}) where {K,T,P}
    empty!(dc.prefix_tree)
    dc.now = zero(T)
    dc.log_likelihood = zero(Float64)
    nothing
end

function copy_clocks!(dst::DirectCall{K,T,P}, src::DirectCall{K,T,P}) where {K,T,P}
    copy!(dst.prefix_tree, src.prefix_tree)
    dst.streams = copy(src.streams)
    return dst
end


"""
    enable!(dc::DirectCall, clock::T, distribution::Exponential, te, when)

Tell the `DirectCall` sampler to enable this clock. The `clock` argument is
an identifier for the clock. The distribution is a univariate distribution
in time. In Julia, distributions are always relative to time `t=0`, but ours
start at some absolute enabling time, ``t_e``, so we provide that here.
The `when` argument is the time at which this clock is enabled, which may be
later than when it was first enabled.

If a particular clock had one rate before an event and it has another rate
after the event, call `enable!` to update the rate.
"""
function enable!(dc::DirectCall{K,T,P}, clock::K, distribution::Exponential,
    te::T, when::T) where {K,T,P}
    dc.prefix_tree[clock] = rate(distribution)
end

function enable!(dc::DirectCall{K,T,P}, clock::K, distribution::D,
    te::T, when::T) where {K,T,P,D<:UnivariateDistribution}
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

function fire!(dc::DirectCall{K,T,P}, clock::K, when::T) where {K,T,P}
    if dc.calculate_likelihood
        dc.log_likelihood += steploglikelihood(dc, dc.now, when, clock)
    end
    disable!(dc, clock, when)
    dc.now = when
end

"""
    next(dc::DirectCall, when::TimeType)

Ask the sampler what clock will be the next to fire and at what time. This does
not change the sampler. You can call this multiple times and get multiple
answers. Each answer is a tuple of `(when, which clock)`. If there is no clock
to fire, then the response will be `(Inf, nothing)`. That's a good sign the
simulation is done.
"""
function next(dc::DirectCall{K,T,P}, when::T) where {K,T,P}
    total = sum!(dc.prefix_tree)
    if total > eps(when)
        rng = race_stream(dc.streams)
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
    return keys(dc.prefix_tree.index)
end

function Base.length(dc::DirectCall)
    return length(dc.prefix_tree)
end

# Implements the interface to return a set of enabled clock keys.
enabled(dc::DirectCall{K,T,P}) where {K,T,P} = enabled(dc.prefix_tree)

isenabled(dc::DirectCall{K,T,P}, clock::K) where {K,T,P} = isenabled(dc.prefix_tree, clock)


function steploglikelihood(dc::DirectCall, now, when, which)
    total = sum!(dc.prefix_tree)
    Δt = when - now
    λ = dc.prefix_tree[which]
    return log(λ) - total * Δt
end

function pathloglikelihood(dc::DirectCall, endtime)
    last_part = if endtime > dc.now
        total = sum!(dc.prefix_tree)
        Δt = endtime - dc.now
        -total * Δt
    else
        zero(Float64)
    end
    return dc.log_likelihood + last_part
end

function Base.haskey(dc::DirectCall{K,T,P}, clock::K) where {K,T,P}
    return isenabled(dc.prefix_tree, clock)
end


Base.haskey(dc::DirectCall{K,T,P}, clock) where {K,T,P} = false
