# Gen.jl and CompetingClocks.jl

From [Gen's](https://www.gen.dev/) point of view, there are two main “hooks” CompetingClocks.jl can plug into:

1. probability distributions (`Distribution{T}` with `random`/`logpdf`), and
2. generative functions (things that implement the Generative Function Interface and produce traces with scores).

CompetingClocks already the operations Gen needs: forward simulation of event paths, and step‑ or path‑wise log‑likelihoods.

Below is an overview of the main integration modes.

---

## 1. CompetingClocks as a Gen distribution over event paths

### What CompetingClocks gives you

Relevant features:

* A `SamplingContext`/`SamplerBuilder` where you can turn on

  * `step_likelihood=true` (for `steploglikelihood`) and/or
  * `path_likelihood=true` (for `pathloglikelihood`), optionally with `likelihood_cnt` for multiple distributions / IS.
* `TrajectoryWatcher` and `PathLikelihoods` objects that do not sample but compute the log‑likelihood of a full path of `(event, time)` pairs. 
* An integration pattern where you replay a sequence of events and call `pathloglikelihood(sampler, end_time)` to get the log probability of that entire trajectory.

This is exactly what Gen’s distribution API wants as `random` (sampling a path) and `logpdf` (log‑likelihood of a path).

### How this looks conceptually in Gen

1. Define a concrete Julia type for a path, e.g. a vector of records `(time, clock_key)` where `clock_key` is whatever you already use.

2. Implement a custom `Distribution{PathType}` (e.g., `GsmpPath`) for Gen:

   * `random(::GsmpPath, params, T, rng)`

     * Build a `SamplingContext` with `path_likelihood=true`.
     * Run your usual integration loop (`initialize_events!`, `next`, `fire!`, `handle_event!`) up to horizon `T`, collecting `(time, key)` into a path.
     * Return that path.
   * `logpdf(::GsmpPath, path, params, T)`

     * Rebuild the same sampler, but now *replay* the given path (possibly using `TrajectoryWatcher`/`PathLikelihoods` instead of a sampler).
     * After replay, return `pathloglikelihood(sampler_or_watcher, T)`, which is a `Float64` when `likelihood_cnt == 1`.

3. In a Gen model you then treat the whole continuous‑time trajectory as one random choice:

   * As a latent path:

     * `path ~ gsmp_path(params, T)` inside an `@gen` function.
   * As an observed path:

     * Use `generate(model, (params, T), choicemap((:path, observed_path)))`; Gen will call your `logpdf` to score that trace.

Once you do this, all Gen inference methods that operate on generative functions calling distributions—importance sampling, SMC, generic MH, etc.—can work with your GSMP path as just another random variable.

---

## 2. CompetingClocks as a custom Generative Function

If you want finer control over the trace (addressable choices per event, incremental updates, gradient hooks), you can wrap CompetingClocks as a full Gen *generative function* rather than just a single distribution.

### Core idea

Gen’s generative function should:

* take simulation parameters and possibly a time horizon as arguments,
* internally run that loop using a `SamplingContext`,
* construct a trace whose choice map encodes the discrete events and/or times,
* set its score to the path log‑likelihood from CompetingClocks.

Key points:

1. **Trace structure and addresses**

   * Gen allows arbitrary address types (not just symbols), and your doc explicitly suggests using a `Pair` key like `:i => 37` “for use with Gen.jl”. 
   * You can align your CompetingClocks clock keys with Gen trace addresses directly:

     * e.g. event key `ClockKey(:infect, i, j)` maps to Gen address `(:infect, i, j)` or to the suggested `:infect => i` / similar scheme.
   * The generative function’s `get_choices(trace)` then gives a choice map whose addresses *are* your event identifiers.

2. **Implementing the Generative Function Interface**

   At minimum you implement:

   * `simulate`: build the `SamplingContext` (with `path_likelihood=true`), run the main loop, record events into the trace, and store `pathloglikelihood` as the trace score.
   * `get_args`, `get_retval` (e.g. the final state or full path), `get_choices`, `get_score`.

   If you want full use of Gen’s MCMC/SMC:

   * add `generate` (to create traces that satisfy constraints on some events/times),
   * `update`/`regenerate` (for efficient incremental proposals that adjust only part of the path),
   * optionally `choice_gradients` and `accumulate_param_gradients!` if you expose derivatives of the path log‑likelihood w.r.t. parameters or event times.

3. **When this is attractive**

   * You want to expose internal events as addressable random choices for custom proposals.
   * You need incremental updates of long trajectories inside MCMC (rather than resimulating whole paths).
   * You want Gen’s gradient‑based machinery to see through to your likelihood (after you supply `logpdf_grad`/`choice_gradients`).

---

## 3. Using CompetingClocks as a likelihood engine inside Gen inference

Even without fully wrapping your simulator as a distribution or generative function, Gen inference code can delegate likelihood calculations to CompetingClocks.

### a) Observation likelihood for event data

Given observed event times/types:

* Represent the data as a list of `(clock_key, time)` events.
* Build a `TrajectoryWatcher` or a sampler created with `path_likelihood=true` and `likelihood_cnt=1`.
* In a Gen model, after sampling parameters, call a pure Julia function that:

  * replays the observed sequence into the watcher/sampler, and
  * returns `log_prob = pathloglikelihood(watcher_or_sampler, end_time)`.

The clean way to plug this into Gen is still via a custom distribution or generative function whose `logpdf` / `get_score` calls this function. That keeps the likelihood inside Gen’s bookkeeping, so all standard inference algorithms see the correct score.

### b) Importance sampling and mixture proposals

We have already seen [importance sampling](../factory.md) with CompetingClocks.

Within Gen:

* A custom inference procedure (e.g. an SMC or MH kernel) can:

  * draw candidate paths using a biased CompetingClocks sampler,
  * use `PathLikelihoods` to evaluate both the target and proposal path log‑likelihoods,
  * compute importance weights or Metropolis–Hastings acceptance ratios, and
  * attach the resulting weight as the trace score (or as an incremental factor inside a custom generative function).

This lets you reuse your variance‑reduction machinery (common random numbers, mixtures, control variates) inside Gen’s programmable inference.

### c) HMC/gradient‑based parameter inference

Two ways this interacts with Gen:

1. **External HMC over latent event times**, Gen for everything else

   * Use your HMC driver to update a latent event list given parameters.
   * Wrap that event list as a variable in a Gen model (e.g. through a distribution over paths).
   * Combine Gen’s inference for other latent variables/parameters with your HMC over the CTDES part.

2. **Gen‑native gradients**

   * Because `pathloglikelihood` is ordinary Julia code, it can in principle be differentiated using AD; Gen’s distribution and generative‑function APIs have hooks (`logpdf_grad`, `choice_gradients`) to expose these gradients to its optimizers and variational methods.

---

## Summary

At a high level, the main interaction patterns are:

1. **Custom Gen distribution over GSMP paths**, with `random` and `logpdf` implemented via CompetingClocks simulation and `pathloglikelihood`.
2. **Custom Gen generative function** that internally runs your simulation, uses your key types as Gen addresses, and sets its score from your path likelihood.
3. **Use of CompetingClocks as a likelihood / weighting engine** in Gen inference code: observation likelihoods for event data, importance‑sampling/mixture proposals (with `PathLikelihoods` and CRNs), and HMC‑style gradient‑based updates driven by `pathloglikelihood`.

All three rely on the same core feature set you already built: configurable samplers/contexts, step‑ and path‑wise log‑likelihoods, watchers, and importance‑sampling tools.

## References

 1. https://www.gen.dev/docs/stable/ref/modeling/distributions/ "Probability Distributions · Gen.jl"
 2. https://www.gen.dev/docs/stable/how_to/custom_gen_fns/ "Adding New Generative Functions · Gen.jl"
 3. https://www.gen.dev/docs/stable/ref/core/gfi/ "Generative Function Interface · Gen.jl"
 4. https://www.gen.dev/docs/dev/how_to/custom_distributions/ "Adding New Distributions · Gen.jl"
