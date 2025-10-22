using Distributions
using Roots
using SpecialFunctions
using Random
using LogExpFunctions: log1mexp

export TranscriptionRate

"""
    expm1x(x)

Compute (exp(x) - 1) / x numerically stably.
For small |x|, uses the Taylor series to avoid catastrophic cancellation.
"""
function expm1x(x::Real)
    if abs(x) < 1e-5
        # Use Taylor series: (exp(x) - 1)/x ≈ 1 + x/2 + x^2/6 + x^3/24 + ...
        return 1.0 + x * (0.5 + x * (1.0 / 6.0 + x * (1.0 / 24.0 + x / 120.0)))
    else
        return expm1(x) / x
    end
end

"""
    _base_cumulative_hazard(α_max, k_rem, t)

Calculates the UN-SHIFTED cumulative hazard Λ(t) = ∫[0, t] λ(u) du.
"""
function _base_cumulative_hazard(α_max::Real, k_rem::Real, t::Real)
    if t < 0
        return 0.0
    end

    # Calculates: α_max * (t - (1 - exp(-k_rem * t)) / k_rem)
    # Uses expm1x for numerical stability near t=0
    integral_term = t * expm1x(-k_rem * t)
    return α_max * (t - integral_term)
end

"""
    _rand_base(rng, α_max, k_rem)

Samples a random time `T` from the UN-SHIFTED (t0=0) distribution
using inverse transform sampling.
"""
function _rand_base(rng::AbstractRNG, α_max::Real, k_rem::Real)
    # We are solving Λ(T) = E
    E = randexp(rng)
    #    f(T) = Λ(T) - E = 0
    f(t) = _base_cumulative_hazard(α_max, k_rem, t) - E

    # 3. Find the root. We know it must be > 0.
    # f(0) = 0 - E = -E < 0, and f(t) → ∞ as t → ∞
    #
    # For a good initial guess, note that as t → ∞:
    # Λ(t) ≈ α_max * (t - 1/k_rem) ≈ α_max * t
    # So t ≈ E / α_max is a reasonable starting point for large E.
    # For small E, we bracket the search between 0 and a reasonable upper bound.
    t_initial_guess = max(E / α_max, 1.0 / k_rem)

    T_base = find_zero(f, t_initial_guess)

    return T_base
end


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
    if t < 0
        return 0.0
    end
    # Forward shift: hazard continues from where it was at t0
    # Λ(t) = ∫[0,t] λ(s) ds = ∫[0,t] α_max * (1 - exp(-k_rem * (s + t0))) ds
    #      = Λ_base(t + t0) - Λ_base(t0)
    return _base_cumulative_hazard(d.α_max, d.k_rem, t + d.t0) -
           _base_cumulative_hazard(d.α_max, d.k_rem, d.t0)
end

function Distributions.cdf(d::TranscriptionRate, t::Real)
    # CDF is F(t) = 1 - exp(-Λ(t))
    if t < 0.0
        return 0.0
    end
    Λ = cumulative_hazard(d, t)
    return 1.0 - exp(-Λ)
end

function Distributions.pdf(d::TranscriptionRate, t::Real)
    # PDF is f(t) = λ(t) * S(t)
    # where λ(t) = α_max * (1 - exp(-k_rem * (t + t0)))
    if t < 0.0
        return 0.0
    end

    # Hazard rate: λ(t) = α_max * (1 - exp(-k_rem * (t + t0)))
    λ = d.α_max * (1.0 - exp(-d.k_rem * (t + d.t0)))

    # Survival function: S(t) = exp(-Λ(t))
    Λ = cumulative_hazard(d, t)
    S = exp(-Λ)

    return λ * S
end

function Distributions.logpdf(d::TranscriptionRate, t::Real)
    if t < 0.0
        return -Inf
    end

    # log(λ(t)) where λ(t) = α_max * (1 - exp(-k_rem * (t + t0)))
    # Use log1mexp for numerical stability
    # log(1 - exp(-x)) = log1mexp(-x) for x > 0
    x = d.k_rem * (t + d.t0)
    log_λ = log(d.α_max) + log1mexp(-x)

    # log(S(t)) = -Λ(t)
    log_S = -cumulative_hazard(d, t)

    return log_λ + log_S
end

"""
    Base.rand(rng::AbstractRNG, d::TranscriptionRate)

Sample from the forward-shifted distribution.
We need to solve: Λ(T) = Λ_base(T + t0) - Λ_base(t0) = E
where E ~ Exp(1).

This is equivalent to: Λ_base(T + t0) = E + Λ_base(t0)
"""
function Base.rand(rng::AbstractRNG, d::TranscriptionRate)
    # Sample E ~ Exp(1)
    E = randexp(rng)

    # We need to solve: Λ_base(T + t0) = E + Λ_base(t0)
    # Let U = T + t0, then we solve: Λ_base(U) = E + Λ_base(t0)
    # This is equivalent to sampling from a base distribution with
    # shifted exponential: E' = E + Λ_base(t0)
    E_shifted = E + _base_cumulative_hazard(d.α_max, d.k_rem, d.t0)

    # Solve Λ_base(U) = E_shifted
    f(u) = _base_cumulative_hazard(d.α_max, d.k_rem, u) - E_shifted

    # Initial guess: start from t0 since U = T + t0 and T ≥ 0
    u_initial_guess = max(E_shifted / d.α_max, d.t0 + 1.0 / d.k_rem)

    U = find_zero(f, u_initial_guess)

    # T = U - t0
    return U - d.t0
end
