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
ordinary `SamplingContext` API. Loading ForwardDiff activates a package
extension; no other setup is required.

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

clock_dist(k, θ) = k == :a ? Exponential(1 / θ[1]) : Exponential(1 / θ[2])

function loglik(θ)
    builder = SamplerBuilder(Symbol, Float64;
        path_likelihood = true,
        likelihood_eltype = eltype(θ),
        method = FirstToFireMethod())
    ctx = SamplingContext(builder, Xoshiro(1))
    enable!(ctx, :a, clock_dist(:a, θ))
    enable!(ctx, :b, clock_dist(:b, θ))
    for (t, k) in trace
        fire!(ctx, k, t)
        enable!(ctx, k, clock_dist(k, θ))
    end
    return pathloglikelihood(ctx, trace[end][1])
end

score = ForwardDiff.gradient(loglik, [0.7, 0.4])
info  = -ForwardDiff.hessian(loglik, [0.7, 0.4])   # observed information
```

For this model the score has a closed form,
``\partial \ell / \partial \lambda_a = n_a/\lambda_a - t_N`` with ``n_a`` the
count of `:a` firings and ``t_N`` the last event time, so the answer is
checkable by hand: `[3/0.7 - 4.5, 2/0.4 - 4.5] ≈ [-0.2143, 0.5]`. Setting the
score to zero gives the maximum-likelihood estimates
``\hat\lambda_a = n_a/t_N`` — the textbook answer, recovered through the
general machinery.

Three things to notice:

- **`likelihood_eltype = eltype(θ)`** is what routes dual numbers into the
  likelihood accumulator. The same closure then works for plain evaluation
  (`eltype` is `Float64`) and differentiation (`eltype` is a dual type),
  including nested duals for `ForwardDiff.hessian`. Setting `likelihood_eltype`
  without requesting a likelihood is an `ArgumentError`.
- **Rebuild the distributions from `θ` inside the function.** The derivative
  information rides in the distribution parameters, so distributions must be
  created from the vector ForwardDiff perturbs. Event times stay `Float64`.
- **`method = FirstToFireMethod()`** (or another general sampler) is needed
  because the automatic selection under `path_likelihood` chooses the
  exponential-only Direct method.

Non-exponential clocks need nothing special — replace `clock_dist` with, say,
`Weibull(θ[1], θ[2])` and the same code differentiates a genuinely semi-Markov
race. Any distribution whose `logpdf`/`logccdf` are differentiable in its
parameters works.

## How it works: the primal boundary

A `SamplingContext` carries both a sampler and likelihood watchers, and under
differentiation they want different numbers. The watchers keep your
dual-parameterized distributions — that is where the gradient accumulates. The
sampler receives a *primal shadow*: the same distribution rebuilt from the
value parts of the parameters, stripped of derivative information. The package
extension performs that stripping at the single boundary through which every
`enable!` reaches the sampler.

This is not a compromise; it is the statistically correct split. Sampling
happens at a concrete parameter point — there is no such thing as drawing a
dual-valued firing time without choosing a derivative convention for the draw
— so the sampler samples at the value of ``\theta`` while the watchers score
with the duals. Everything the sampler does (heaps of pending times,
common random numbers, cloning) stays in plain `Float64`.

One consequence to know about: with `step_likelihood = true` and a dual
`likelihood_eltype`, the context always installs a tracking watcher rather
than using a sampler's native step-likelihood, because the native computation
would run on primal distributions and return a derivative-free `Float64`.
You do not need to do anything; the builder handles it.

## Walk up to an event, then sample — with the score in hand

Because the sampler stays functional, the context's walk-up promise survives
differentiation: replay an observed prefix, then let the sampler continue.

```julia
function continuation(θ)
    builder = SamplerBuilder(Symbol, Float64;
        path_likelihood = true,
        likelihood_eltype = eltype(θ),
        method = FirstToFireMethod())
    ctx = SamplingContext(builder, Xoshiro(7))
    enable!(ctx, :a, clock_dist(:a, θ))
    enable!(ctx, :b, clock_dist(:b, θ))
    for (t, k) in trace[1:3]           # replay the observed prefix
        fire!(ctx, k, t)
        enable!(ctx, k, clock_dist(k, θ))
    end
    (when, which) = next(ctx)          # sample a continuation at primal θ
    fire!(ctx, which, when)
    return pathloglikelihood(ctx, when)
end

ForwardDiff.gradient(continuation, [0.7, 0.4])
```

The continuation is sampled from the primal measure while the log-likelihood
stays dual, which is exactly the likelihood-ratio construction: for a
functional ``f`` of the sampled continuation,
``f \cdot \nabla_\theta \log L`` estimates ``\nabla_\theta E[f]`` without
bias. The context is not just a trace scorer; it is a score-function gradient
estimator.

## Differentiating the next-event probabilities

`stepconditionalprobability(tw, t)` maps out ``P[K \mid T = t]``, the
probability that each enabled clock is the one that fires at time ``t``.
Its entries are hazard ratios, and they differentiate the same way — enable
dual-parameterized distributions on a `TrackWatcher` and read out a chosen
clock's probability. (Both names live on the unexported developer surface,
so qualify them.)

```julia
function p_a_fires(θ)
    tw = CompetingClocks.TrackWatcher{Symbol,Float64}()
    enable!(tw, :a, Weibull(θ[1], θ[2]), 0.0, 0.0, Xoshiro(1))
    enable!(tw, :b, Exponential(1 / θ[3]), 0.0, 0.0, Xoshiro(1))
    return CompetingClocks.stepconditionalprobability(tw, 1.3)[:a]
end

ForwardDiff.gradient(p_a_fires, [1.6, 2.0, 0.4])
```

Because the probabilities sum to one, the gradients across clocks sum to
zero — a useful check on any model where you use these sensitivities.

## The low-level route

The watchers work standalone if you want no context at all: construct
`TrajectoryWatcher{K,T,L}()` with `L = eltype(θ)` and drive
`enable!`/`fire!`/`pathloglikelihood` directly, or hold your own accumulator
and sum [`steploglikelihood`](@ref) over a `TrackWatcher` for per-step scores.
The two-parameter constructors `TrajectoryWatcher{K,T}()` and
`PathLikelihoods{K,T}(cnt)` remain `Float64`. `PathLikelihoods{K,T,L}` scores
a whole vector of candidate parameterizations in one replay — useful when
tuning importance-sampling mixtures (see
[Importance Sampling](importance_skills.md)).

## What this does and does not differentiate

This page differentiates the likelihood of a *fixed* trace (plus, in the
walk-up form, sampled continuations weighted by their scores). That is exactly
what parameter inference wants: the data pin the event order, and the
likelihood is smooth in ``\theta``. It is **not** a pathwise derivative
*through the simulation*: if you change ``\theta`` and re-run the sampler,
which clock wins each race can change discontinuously, and the expectation of
a simulation output needs estimators that account for those event-order
changes (pathwise/IPA where valid, the score-function estimator built from
exactly the gradient computed here, or weak-derivative methods in general).

For using path likelihoods with Hamiltonian Monte Carlo over event times, see
[Hamiltonian Monte Carlo](hamiltonianmontecarlo.md).
