# Integration with Survival.jl

Fitting a survival distribution and using it as a GSMP clock.

This example shows how to:

1. Generate some noisy event times with random right‐censoring.
2. Fit a nonparametric Kaplan–Meier survival curve using `Survival.jl`.
3. Fit a parametric Weibull distribution to that curve.
4. Wrap the fitted Weibull as a `UnivariateDistribution` and use it as a
   clock distribution in `CompetingClocks.jl`.

```julia
using Random
using Distributions
using Survival          # Kaplan–Meier and EventTime
using LsqFit            # nonlinear least squares for parametric fit
using CompetingClocks   # GSMP sampler

Random.seed!(2024)

# ------------------------------------------------------------
# 1. Synthetic "messy" event times with right censoring
# ------------------------------------------------------------

n = 200
true_dist = Weibull(1.7, 4.0)             # ground-truth failure-time distribution

t_event  = rand(true_dist, n)             # actual failure times
t_censor = rand(Uniform(0.0, 8.0), n)     # censoring times

time   = min.(t_event, t_censor)          # observed follow-up time
status = t_event .<= t_censor             # true = failure observed, false = censored

# Optionally store as EventTime values (useful elsewhere in Survival.jl)
event_times = EventTime.(time, status)

# ------------------------------------------------------------
# 2. Nonparametric Kaplan–Meier survival estimate
# ------------------------------------------------------------

km = fit(KaplanMeier, time, status)

t_km = km.events.time     # unique event/censor times
S_km = km.survival        # Kaplan–Meier Ŝ(t)

# ------------------------------------------------------------
# 3. Parametric Weibull fit to the Kaplan–Meier curve
#    S_model(t; α, θ) = P(T > t) = ccdf(Weibull(α, θ), t)
# ------------------------------------------------------------

function km_model(t, p)
    d = Weibull(p[1], p[2])              # p[1] = α (shape), p[2] = θ (scale)
    return [ccdf(d, ti) for ti in t]     # model survival at each time ti
end

p0 = [1.0, 3.0]                          # initial guess for (α, θ)
fit_result = curve_fit(km_model, t_km, S_km, p0)

α̂, θ̂ = fit_result.param                # fitted Weibull parameters
failure_time_dist = Weibull(α̂, θ̂)      # this is a ContinuousUnivariateDistribution

# ------------------------------------------------------------
# 4. Use the fitted distribution as a CompetingClocks clock
# ------------------------------------------------------------

const ClockKey = Symbol

builder = SamplerBuilder(ClockKey, Float64)
rng = Xoshiro(1234)
sampler = SamplingContext(builder, rng)

# Enable a GSMP clock whose waiting time is given by the fitted Weibull
enable!(sampler, :failure, failure_time_dist)

# Sample the next event time from this calibrated clock
when, which = next(sampler)
@show when which
```
