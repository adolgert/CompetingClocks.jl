using Distributions: UnivariateDistribution

"""
    Delayed{D1,D2}(initiation, duration)

Pair of distributions for delayed reactions.

- `initiation` : distribution for time from enabling to initiation
- `duration`   : distribution for time from initiation to completion

At the package's API boundary you may pass a plain `Pair` of distributions,
`initiation => duration`, which is converted to a `Delayed`. For example
`enable!(ctx, clock, Exponential(1.0) => Gamma(2, 0.5))`.

Note: `Delayed` is not itself a distribution. It's a marker type that triggers
specialized `enable!` behavior in `SamplingContext`.
"""
struct Delayed{D1<:UnivariateDistribution,D2<:UnivariateDistribution}
    initiation::D1
    duration::D2
end

"""
    Delayed(p::Pair{<:UnivariateDistribution,<:UnivariateDistribution})

Construct a `Delayed` from an ordinary `Pair` of distributions,
`initiation => duration`. This lets callers write the pair syntax
`Exponential(1.0) => Gamma(2, 0.5)` without any type piracy on `Base.:(=>)`.
"""
Delayed(p::Pair{<:UnivariateDistribution,<:UnivariateDistribution}) = Delayed(p.first, p.second)


"""
    DelayedState{K,T}

Holds per-clock delayed state for a `SamplingContext` whose user key type is `K`
and time type is `T`.

`durations` maps a user clock key to the **duration distribution** for the
completion phase.
"""
mutable struct DelayedState{K,T}
    durations::Dict{K,UnivariateDistribution}
end

"""
    DelayedState{K,T}()

Construct an empty delayed state.
"""
DelayedState{K,T}() where {K,T} = DelayedState{K,T}(Dict{K,UnivariateDistribution}())


"""
    reset!(ds::DelayedState)

Clear all stored duration distributions.
"""
function reset!(ds::DelayedState)
    empty!(ds.durations)
    return ds
end

"""
    clone(ds::DelayedState)

Deep copy of a `DelayedState`.
"""
function clone(ds::DelayedState{K,T}) where {K,T}
    return DelayedState{K,T}(copy(ds.durations))
end
