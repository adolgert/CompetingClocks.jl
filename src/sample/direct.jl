using DataStructures
using Random: rand, AbstractRNG
using Distributions: Uniform, Exponential, rate

export DirectCall, enable!, disable!, next


"""
    DirectCall{T}

DirectCall is responsible for sampling among Exponential distributions. It
samples using the Direct method. In this case, there is no optimization to
that Direct method, so we call it DirectCall because it recalculates
everything every time you call it.

The type `T` is the type of an identifier for each transition. This identifier
is usually a nominal integer but can be a any key that identifies it, such as
a string or tuple of integers. Instances of type `T` are used as keys in a
dictionary.
"""
struct DirectCall{T,P}
    # Map from clock name to index in propensity array.
    index::Dict{T, Int64}
    # Map from index in propensity array to clock name.
    key::Vector{T}
    prefix_tree::P
    DirectCall{T,P}(tree::P) where {T,P} = new(
        Dict{T, Int64}(), Vector{T}(), tree
        )
end

DirectCall{T}() where {T} =
    DirectCall{T,CumSumPrefixSearch{Float64}}(CumSumPrefixSearch(Float64))


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
function enable!(dc::DirectCall{T}, clock::T, distribution::Exponential,
        te::Float64, when::Float64, rng::AbstractRNG) where {T}
    hazard = rate(distribution)
    idx = get(dc.index, clock, 0)
    if idx == 0
        push!(dc.prefix_tree, hazard)
        dc.index[clock] = length(dc.prefix_tree)
        push!(dc.key, clock)
    else
        dc.prefix_tree[idx] = hazard
    end
end


"""
    disable!(dc::DirectCall, clock::T, when)

Tell the `DirectCall` sampler to disable this clock. The `clock` argument is
an identifier for the clock. The `when` argument is the time at which this
clock is enabled.
"""
function disable!(dc::DirectCall{T}, clock::T, when::Float64) where {T}
    dc.prefix_tree[dc.index[clock]] = 0.0
end


"""
    next(dc::DirectCall, when::Float64, rng::AbstractRNG)

Ask the sampler what clock will be the next to fire and at what time. This does
not change the sampler. You can call this multiple times and get multiple
answers. Each answer is a tuple of `(when, which clock)`. If there is no clock
to fire, then the response will be `(Inf, nothing)`. That's a good sign the
simulation is done.
"""
function next(dc::DirectCall, when::Float64, rng::AbstractRNG)
    if length(dc.prefix_tree) == 0
        return (Inf, nothing)
    end
    total = sum!(dc.prefix_tree)
    if total > eps(Float64)
        chosen, hazard_value = rand(rng, dc.prefix_tree)
        @assert chosen < length(dc.prefix_tree) + 1
        tau = when + rand(rng, Exponential(1 / total))
        return (tau, dc.key[chosen])
    else
        return (Inf, nothing)
    end
end
