# Changelog

## v0.3.1 (unreleased)

### Added

- **The path log-likelihood is differentiable with ForwardDiff.** Once a trace
  is fixed, `pathloglikelihood` is a smooth function of the clock
  distributions' parameters, and its gradient is the score. The user-facing
  entry point is a new `SamplerBuilder` keyword, `likelihood_eltype`: build the
  context with `likelihood_eltype=eltype(θ)` inside a closure and
  `ForwardDiff.gradient`/`hessian` work through `pathloglikelihood` and
  `steploglikelihood`, including nested duals for observed information. Setting
  it without requesting a likelihood is an `ArgumentError`. See the new manual
  page "Differentiating the Likelihood".
- **Primal boundary between watchers and samplers.** The likelihood watchers
  keep the caller's (possibly dual-parameterized) distributions while the
  sampler receives a value-only primal copy through
  `primal_distribution` — the identity for ordinary parameters, with a
  ForwardDiff package extension (`[weakdeps]`, loaded only when ForwardDiff
  is) that strips dual parameters, recursing for nested duals and
  reconstructing `truncated` wrappers. Consequently `next(ctx)` keeps working
  under differentiation: it samples a continuation at the primal parameter
  point while the accumulated log-likelihood stays dual — replaying an
  observed prefix, sampling forward, and reading the dual `pathloglikelihood`
  is the score-function estimator `∇E[f] = E[f·∇log L]`. With
  `step_likelihood=true` and a non-`Float64` eltype the context installs a
  `TrackWatcher` even when the sampler computes step likelihoods natively,
  because the native path would run on primal distributions and silently
  return an underived `Float64`.
- `TrajectoryWatcher{K,T,L}` and `PathLikelihoods{K,T,L}` gained the
  accumulator number type `L` as a third type parameter; the two-parameter
  constructors remain and mean `L=Float64`.
- `hazard` accepts `Real` times, and `stepconditionalprobability` — the
  next-event probabilities `P[K | T=t]` — is differentiable: its `Dict` value
  type follows the distributions' parameter number type (`Float64` for plain
  distributions, unchanged). Its `t < te` early return promotes to the
  parameter type so a dual query never evaluates `logpdf` at a negative age.

### Changed

- On the low-level (unexported, `public`) surface, the watcher type identity
  changed: `typeof(tw) == TrajectoryWatcher{K,T}` comparisons are now false
  (the concrete type is `TrajectoryWatcher{K,T,Float64}`), watcher objects
  serialized under the two-parameter layout do not deserialize, and a concrete
  field annotation `::TrajectoryWatcher{K,T}` is now abstract. Dispatch
  constraints and constructors are unaffected. `SamplerBuilder`'s struct
  gained a field, so positional construction (not the keyword API) changed
  arity.

### Internal

- ForwardDiff is a test-only and docs-only dependency plus a weak dependency
  for the extension; the package's hard dependencies are unchanged.

## v0.3.0

### Breaking

- **`fire!` and `disable!` are now different operations.** Firing realizes a
  clock's draw (fully consumed; a re-enabled clock draws fresh); disabling
  censors it (a sampler may retain the remaining randomness for reuse).
  `CombinedNextReaction.fire!` now consumes the draw, and its `next` is pure —
  the survival mark that previously ran inside `next` moved to `fire!`.
  `MultiSampler.fire!` forwards `fire!` (not `disable!`) to the owning
  sub-sampler. This is breaking in the semantic sense for low-level users who
  relied on `fire! ≡ disable!` on `CombinedNextReaction`, or who called
  `next()` expecting it to mark the returned clock. In-contract high-level
  users see identical sample paths. Low-level users must call `fire!` for the
  clock that fired and `disable!` only for clocks whose preconditions
  vanished; see the manual's "Re-enabling and Memory" page.
- **The exported API is the high-level layer.** `SamplingContext`, the
  `SamplerBuilder`/spec machinery, `Delayed`, and `Never` are exported. The
  low-level layer (concrete samplers, watchers, `SSA`, `MultiSampler`, common
  random number internals) remains fully supported for framework authors via
  qualified `CompetingClocks.X` access and is declared `public` on Julia 1.11
  and later; it is no longer exported. `getindex`, `keys`, and `length` are no
  longer exported (the `Base` methods dispatch without the binding).
- **`Base.:(=>)` is no longer overloaded on distributions** (it was type
  piracy that changed `Pair` semantics for unrelated code). The user-facing
  syntax is unchanged: `enable!(ctx, clock, initiation => duration)` now works
  by dispatch on an ordinary `Pair`; `Delayed(d1, d2)` remains the explicit
  constructor. A `Pair` on a context built without `support_delayed=true`
  errors at the boundary.
- On a delayed context, plain `next(ctx)` errors with a pointer to
  `next_delayed`, and the tuple `(clock, phase)` is documented as the public
  event identity.

### Added

- Documented contract for `next`: the return value is a *reservation*, valid
  until the next `enable!`/`disable!`/`fire!`; the time argument never
  decreases and never advances past a pending firing without that event being
  fired (the last fired event's time is always safe). `MultiSampler`'s
  superposition satisfies the contract by construction.
- `fire!` is declared on the abstract `SSA` interface with a documented
  fallback to `disable!` — correct for every sampler that retains no residual
  draw randomness; override it (as `CombinedNextReaction` does) when firing
  must consume retained randomness.
- Statistical correctness gauntlet under `test/gauntlet/`: Doob–Meyer
  step-cumulant, mark calibration, and a contract-conforming two-sample
  Anderson–Darling comparison against `FirstReaction` as the trusted
  reference, run over a sampler-by-condition matrix (exponential, Weibull,
  and shifted-enabling conditions; forget/remember memory; cycle/complete
  graphs) with Benjamini–Hochberg control across the matrix. The current
  matrix shows no statistically supported defect in any sampler.
- Manual: a Concepts section (background, notation, GSMPs, shifted sampling,
  vector addition systems) and a new "Re-enabling and Memory" page covering
  the three re-enabling idioms and the fire-versus-disable distinction.
- Crosswalk table between sampler specs (`FirstToFireMethod`, ...) and the
  underlying sampler types.

### Fixed

- `PSSACR`: `copy_clocks!` and the cached-`next` path read a nonexistent
  field of `OrderedSample`; copying with a populated cache threw.
- `MultiSampler.next` queried the second and later sub-samplers with the
  previous sub-sampler's firing time instead of the true current time.
- `MultiSampler.fire!` routed to the sub-sampler's `disable!`, which skipped
  `DirectCall`'s native likelihood accumulation.
- `SamplingContext` capability probes used `catch MethodError` (which binds
  any exception rather than filtering); they now use the
  `has_steploglikelihood`/`has_pathloglikelihood` traits, and errors from
  sampler likelihood code propagate. The `has_pathloglikelihood` trait
  definitions missing `<:` never matched concrete types.
- `MemorySampler` now forwards `fire!` as its docstring promised.
- Removed the exported-but-undefined `ChatReaction`; fixed the
  `DirectMethod` docstring advertising a `:scan` option the constructor
  rejects; standardized `SamplerChoice{SamplerKey,Key}` parameter order
  between `MultiSampler` and `MultipleDirect`.
- The `enable!` docstring's `truncated()` guidance now pairs truncation with
  the shifted enabling time; truncation alone displaces every firing by the
  truncation point.
- Removed unused dependencies (`Documenter`, `Combinatorics`,
  `SpecialFunctions`, `InteractiveUtils`); `available_samplers()` no longer
  requires the REPL stdlib to be loaded.

### Internal

- `SamplingContext` middleware (likelihood, debug/recording) is carried in a
  single compile-time watcher tuple; allocations and inference are unchanged,
  and adding a watcher no longer touches every `enable!`/`fire!` overload.

## v0.2.0

- Simplified interface to samplers; path likelihoods; delayed reactions;
  Rejection-based SSA and partial-propensity composition-rejection samplers;
  importance sampling with mixtures; examples for Gen.jl, Turing.jl, and
  Survival.jl integration; tracing and debugging support. Registered in the
  Julia General registry.
