# Differentiating the Path Likelihood

The log-likelihood of a recorded trajectory is a smooth function of the clock
distributions' parameters. Once the sequence of events is fixed — because you
observed it, or because you recorded it from a simulation run — the path
log-likelihood is a sum of `logpdf` and `logccdf` terms, and its gradient with
respect to the distribution parameters is the *score function*,

```math
s(\theta) = \nabla_\theta \log L(\theta; \text{trace}).
```

The score is the working part of most gradient-based statistics: maximum
likelihood by Newton or quasi-Newton steps, Hamiltonian Monte Carlo over
parameters, observed Fisher information for standard errors, and
likelihood-ratio (score-function) sensitivity estimates of the form
``\nabla_\theta E[f] = E[f \cdot s(\theta)]``.

This page shows how to compute the score with
[ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) through the
watchers in this package. No special mode is required; the likelihood
accumulators are generic over the number type.

## How it works

ForwardDiff evaluates your function with dual numbers — a value carrying
derivative components — in place of `Float64`. Three facts let a dual flow from
``\theta`` to the returned log-likelihood:

1. A distribution built from dual parameters, such as `Weibull(k, s)` with
   dual `k` and `s`, is still a `UnivariateDistribution`, so it passes through
   `enable!` like any other distribution.
2. `logpdf` and `logccdf` of that distribution at a `Float64` time return a
   dual carrying ``\partial \log f / \partial \theta``.
3. The watcher's accumulator type is a type parameter: `TrajectoryWatcher{K,T,L}`
   stores its running log-likelihood as an `L`. Construct it with
   `L = eltype(θ)` and the same code runs in plain `Float64` mode and in dual
   mode. The two-parameter form `TrajectoryWatcher{K,T}()` remains `Float64`.

Event *times* stay `Float64` throughout — you are differentiating with respect
to the parameters, not the times, and a fixed trace pins the event order, which
is what makes ``\theta \mapsto \log L`` smooth.

## Example: the score of a two-clock race

Two clocks race repeatedly; the trace records which fired when. Clock `:a`
and clock `:b` are exponential with rates ``\lambda_a`` and ``\lambda_b``
(note that `Distributions.Exponential` takes the *scale*, the inverse rate).
Each fired clock is re-enabled fresh at its firing time.

```julia
using CompetingClocks
using Distributions
using ForwardDiff
using Random: Xoshiro

# The observed data: (firing time, which clock). Think of a maintenance log.
trace = [(0.9, :a), (1.7, :b), (2.2, :a), (3.1, :a), (4.5, :b)]
rng = Xoshiro(1)  # the interface requires an RNG; the watcher never draws from it

function loglik(θ)
    λa, λb = θ[1], θ[2]
    tw = CompetingClocks.TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
    enable!(tw, :a, Exponential(1 / λa), 0.0, 0.0, rng)
    enable!(tw, :b, Exponential(1 / λb), 0.0, 0.0, rng)
    for (t, k) in trace
        fire!(tw, k, t)
        d = (k == :a) ? Exponential(1 / λa) : Exponential(1 / λb)
        enable!(tw, k, d, t, t, rng)
    end
    return pathloglikelihood(tw, trace[end][1])
end

score = ForwardDiff.gradient(loglik, [0.7, 0.4])
```

For this model the score has a closed form,
``\partial \ell / \partial \lambda_a = n_a/\lambda_a - t_N`` with ``n_a`` the
count of `:a` firings and ``t_N`` the last event time, so the answer is
checkable by hand: `[3/0.7 - 4.5, 2/0.4 - 4.5] ≈ [-0.2143, 0.5]`. ForwardDiff
reproduces it to machine precision. Setting the score to zero gives the
maximum-likelihood estimates ``\hat\lambda_a = n_a/t_N`` — the textbook answer,
recovered through the general machinery.

The pattern has three rules:

- **Construct the watcher with `eltype(θ)`.** That is what routes dual numbers
  into the accumulator.
- **Rebuild the distributions from `θ` inside the function.** The derivative
  information rides in the distribution parameters, so the distributions must
  be created from the vector that ForwardDiff perturbs.
- **Leave times alone.** `te` and `when` arguments remain `Float64`.

## Non-exponential clocks

Nothing above is special to the exponential. A Weibull clock racing an
exponential one — a genuinely semi-Markov race — differentiates the same way:

```julia
function loglik_weibull(θ)
    shape, scale, λb = θ[1], θ[2], θ[3]
    tw = CompetingClocks.TrajectoryWatcher{Symbol,Float64,eltype(θ)}()
    enable!(tw, :a, Weibull(shape, scale), 0.0, 0.0, rng)
    enable!(tw, :b, Exponential(1 / λb), 0.0, 0.0, rng)
    for (t, k) in trace
        fire!(tw, k, t)
        d = (k == :a) ? Weibull(shape, scale) : Exponential(1 / λb)
        enable!(tw, k, d, t, t, rng)
    end
    return pathloglikelihood(tw, trace[end][1])
end

ForwardDiff.gradient(loglik_weibull, [1.6, 2.0, 0.4])
```

Any distribution whose `logpdf`/`logccdf` are differentiable in its parameters
works, which covers the standard parametric families. The left-shift
conditioning terms for clocks enabled with past zero-points (`te < when`) are
part of the same accumulated sum and differentiate with everything else.

## The per-step route

If you want per-step scores — for stochastic-gradient methods, or to weight
individual steps — skip the accumulating watcher and sum
[`steploglikelihood`](@ref) yourself over a `TrackWatcher`, holding your own
accumulator:

```julia
function loglik_steps(θ)
    λa, λb = θ[1], θ[2]
    tw = CompetingClocks.TrackWatcher{Symbol,Float64}()
    enable!(tw, :a, Exponential(1 / λa), 0.0, 0.0, rng)
    enable!(tw, :b, Exponential(1 / λb), 0.0, 0.0, rng)
    total = zero(eltype(θ))
    t0 = 0.0
    for (t, k) in trace
        total += steploglikelihood(tw, t0, t, k)
        fire!(tw, k, t)
        d = (k == :a) ? Exponential(1 / λa) : Exponential(1 / λb)
        enable!(tw, k, d, t, t, rng)
        t0 = t
    end
    return total
end
```

The two routes agree; the per-step form needs no type parameter at all because
the caller owns the accumulator.

## Scoring one trace under many parameterizations

`PathLikelihoods{K,T,L}` is the vector version of the same idea: each
`enable!` may take a vector of candidate distributions, and one pass over the
trace returns one log-likelihood per candidate. Its accumulator takes the same
third type parameter, so a dual-typed `PathLikelihoods` differentiates a whole
family of parameterizations in a single replay — useful when tuning
importance-sampling mixtures (see [Importance Sampling](importance_skills.md)).

## What this does and does not differentiate

This page differentiates the likelihood of a *fixed* trace. That is exactly
what parameter inference wants: the data pin the event order, and the
likelihood is smooth in ``\theta``. It is **not** a derivative *through the
simulation*: if you change ``\theta`` and re-run the sampler, which clock wins
each race can change discontinuously, and the expectation of a simulation
output needs estimators that account for those event-order changes
(pathwise/IPA where valid, score-function using exactly the gradient computed
here, or weak-derivative methods in general). The score you compute on
recorded traces is the unbiased likelihood-ratio building block:
``\nabla_\theta E[f] = E[f \cdot \nabla_\theta \log L]``.

For using path likelihoods with Hamiltonian Monte Carlo over event times, see
[Hamiltonian Monte Carlo](hamiltonianmontecarlo.md).
