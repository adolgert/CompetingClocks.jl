using Distributions
using Roots
using SpecialFunctions
using Random
using LogExpFunctions: log1mexp

export TranscriptionRate


"""
    TranscriptionRate(α_max, k_rem; t0=0.0)

A continuous univariate distribution defined by a time-varying
hazard rate that starts at time `t0` where `t0 > 0`.

The hazard rate λ_s(t) is:
- 0                      if t < 0
- α_max * (1 - exp(-k_rem * (t + t0))) if t >= 0

Parameters:
- `α_max`: The maximum asymptotic hazard rate.
- `k_rem`: The rate parameter controlling how fast the hazard approaches α_max.
- `t0`: The shift time. The distribution's support is [0, ∞).
"""
struct TranscriptionRate{T<:Real} <: ContinuousUnivariateDistribution
    α_max::T
    k_rem::T
    t0::T

    function TranscriptionRate(α_max::Real, k_rem::Real; t0::Real=0.0)
        if α_max <= 0 || k_rem <= 0
            error("Parameters α_max and k_rem must be positive.")
        end
        if t0 < 0
            error("Parameter t0 must be non-negative.")
        end
        T = promote_type(typeof(α_max), typeof(k_rem), typeof(t0))
        new{T}(α_max, k_rem, t0)
    end
end

Distributions.params(d::TranscriptionRate) = (d.α_max, d.k_rem, d.t0)


"""
    cumulative_hazard(d::TranscriptionRate, t::Real)

Calculates the cumulative hazard Λ(t) for the forward-shifted distribution.
The hazard λ(t) = α_max * (1 - exp(-k_rem * (t + t0))) integrates to:
Λ(t) = Λ_base(t + t0) - Λ_base(t0)
"""
function cumulative_hazard(d::TranscriptionRate, t::Real)
    t < zero(t) && return zero(t)
    return d.α_max * (t + exp(-d.k_rem * d.t0) * expm1(-d.k_rem * t) / d.k_rem)
end

function hazard(d::TranscriptionRate, t::Real)
    return d.α_max * (1.0 - exp(-d.k_rem * (t + d.t0)))
end

function Distributions.cdf(d::TranscriptionRate, t::Real)
    t < 0.0 && return 0.0
    return 1.0 - exp(-cumulative_hazard(d, t))
end

function Distributions.pdf(d::TranscriptionRate, t::Real)
    t < 0.0 && return 0.0
    # PDF is f(t) = λ(t) * S(t)
    # where λ(t) = α_max * (1 - exp(-k_rem * (t + t0)))
    survival = exp(-cumulative_hazard(d, t))
    return hazard(d, t) * survival
end

function Distributions.logpdf(d::TranscriptionRate, t::Real)
    t < 0.0 && return -Inf

    # log(λ(t)) where λ(t) = α_max * (1 - exp(-k_rem * (t + t0)))
    # Use log1mexp for numerical stability
    # log(1 - exp(-x)) = log1mexp(-x) for x > 0
    x = d.k_rem * (t + d.t0)
    log_λ = log(d.α_max) + log1mexp(-x)

    # log(S(t)) = -Λ(t)
    log_S = -cumulative_hazard(d, t)

    return log_λ + log_S
end

function Distributions.logccdf(d::TranscriptionRate, t::Real)
    t < 0.0 && return 0.0
    return -cumulative_hazard(d, t)
end

"""
    Base.rand(rng::AbstractRNG, d::TranscriptionRate)

Sample from the forward-shifted distribution.
We need to solve: Λ(T) = Λ_base(T + t0) - Λ_base(t0) = E
where E ~ Exp(1).

This is equivalent to: Λ_base(T + t0) = E + Λ_base(t0)
"""
function Base.rand(rng::AbstractRNG, d::TranscriptionRate)
    E = randexp(rng)
    f(t) = cumulative_hazard(d, t) - E
    if d.t0 < 1e-5
        t_initial_guess = sqrt(2.0 * E / (d.α_max * d.k_rem))
    else
        t_initial_guess = E / (d.α_max * (1.0 - exp(-d.k_rem * d.t0)))
    end
    return find_zero(f, t_initial_guess)
end
