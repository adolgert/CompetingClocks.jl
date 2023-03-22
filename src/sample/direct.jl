using DataStructures
using Random: rand, AbstractRNG
using Distributions: Uniform, Exponential, params

export DirectCall, enable!, disable!, next


"""
    DirectCall{T}

DirectCall is responsible for sampling among Exponential distributions. It
samples using the Direct method. In this case, there is no optimization to
that Direct method.

The type `T` is the type of an identifier for each transition. This identifier
is usually a nominal integer but can be a any key that identifies it, such as
a string or tuple of integers. Instances of type `T` are used as keys in a
dictionary.
"""
struct DirectCall{T}
    key::Dict{T, Int64}
    propensity::Vector{Float64}
    DirectCall{T}() where {T} = new(Dict{T, Int64}(), zeros(Float64, 0))
end


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
    hazard = params(distribution)[1]
    if !haskey(dc.key, clock)
        dc.key[clock] = length(push!(dc.propensity, hazard))
    else
        dc.propensity[dc.key[clock]] = hazard
    end
end


"""
    disable!(dc::DirectCall, clock::T, when)

Tell the `DirectCall` sampler to disable this clock. The `clock` argument is
an identifier for the clock. The `when` argument is the time at which this
clock is enabled.
"""
function disable!(dc::DirectCall{T}, clock::T, when::Float64) where {T}
    dc.propensity[dc.key[clock]] = 0.0
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
    total = sum(dc.propensity)
    if total > eps(Float64)
        chosen = searchsortedfirst(cumsum(dc.propensity), rand(rng, Uniform(0, total)))
        @assert chosen < length(dc.propensity) + 1
        key_name = [x for (x, y) in pairs(dc.key) if y == chosen][1]
        return (when - log(rand(rng)) / total, key_name)
    else
        return (Inf, nothing)
    end
end
