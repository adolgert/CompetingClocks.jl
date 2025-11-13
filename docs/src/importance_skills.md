# Importance Sampling for Simulation

## The Process

When you apply importance sampling in simulation, the workflow feels like this:

1. **Baseline run:** Simulate under the desired distributions, $p(x)$. If the event is rare, it might not ever happen in your simulations or might happen too few times to get good statistics on it.
2. **Bias intuition:** Identify which physical parameters or transitions make the event unlikely (e.g., promoter turns off too soon). Modify those rates to define $q(x)$. Try to pick rates that happen once or twice in the simulation, and start with gentle bias.
3. **Run biased simulations:** Under $q$, the rare event occurs more often.
4. **Compute path likelihoods:** Using both the *biased* and *true* rates, form $w = p(x)/q(x)$ for each path.
5. **Reweight results:** Estimate probabilities or expectations with the weighted average.
6. **Diagnose:** Check the mean and variance of weights, or compute ESS. If weights vary over many orders of magnitude, adjust $q$ to reduce the gap.
7. **Iterate:** Gradually refine the proposal until ESS stabilizes and probability estimates converge.

The main problem is that too large of bias on distributions can lead to mathematical underflow in calculation of the weights. Intuitively, a stochastic simulation can have a lot of individual sampled events, and each event's probability multiplies to get the probability of a path of samples in a trajectory. If those samples are repeatedly biased, they can cause numbers that are too small to represent.
```math
w = \frac{L(\lambda_{\mbox{target}})}{L(\lambda_{\mbox{proposal}})} = \left(\frac{\lambda_{\mbox{target}}}{\lambda_{\mbox{proposal}}}\right)^N e^{-(\lambda_{\mbox{target}} - \lambda_{\mbox{proposal}})T}
```
What you'll see in practice is that the initial simulation, under $p$, works fine, that a small change in a distribution's parameters still works fine, and then the importance-weighted estimates fall off a cliff and show values like $10^{-73}$.

If distributions are scaled, the mean weight can be much less than one. A mean weight that isn't small is a good sign. We want a sense of the ratio of $w=p/q$ the ratio of weights.

## Evaluate proposal quality

### Track how ESS scales with N

```julia
# Stable log-sum-exp
function logsumexp(x)
    m = maximum(x)
    return m + log(sum(exp.(x .- m)))
end

# Estimate E_q[w] = 1 using the *unshifted* Δ (i.e., don't subtract max here)
# log mean weight = log( (1/N) * sum(exp(Δ)) )
log_mean_w = logsumexp(Δ) - log(length(Δ))
mean_w_est = exp(log_mean_w)
println("E_q[w] (should be ~1): ", mean_w_est)

# For self-normalized IS you don’t need this, but it’s a nice sanity check.

println("mean weight ≈ ", mean_w_est)
w = exp.(Δ)                  # unshifted
# normalize stably:
Z_log = logsumexp(Δ)
wn = exp.(Δ .- Z_log)        # now sum(wn) = 1
ESS_norm = 1 / sum(wn.^2)    # between 1 and N
# Close to (sum(importance)^2) / sum(importance.^2)
println("ESS (normalized): ", ESS_norm)
```
ESS is explained sum of squares. You're looking for this to be a significant fraction of $N$ the number of trajectories.

If it stops growing linearly with N, your proposal is too far from the target.

### Use the coefficient of variation of weights to see where variance explodes.

```math
\mbox{CV}^2 = \frac{\mbox{Var}(w)}{E[w]^2} = \frac{N}{\mbox{ESS}} - 1
```
Coefficient of variation is a reparametrization of ESS.

### Make a weight histogram

Plot log of the weight, $\log_{10}(w)$ or $\log(w/\mbox{mean}(w))$. You want a unimodal, not-too-wide shape. Heavy-tailed distributions indicate you're close to degeneracy.

## Proposal Improvement

### Mixture proposals

Instead of one biased model, use a mixture of proposals. This helps if different regions of state space are rare for different reasons. For the gene example, maybe one proposal distribution keeps the promoter on longer and another emphasizes instead reducing degradation of MRNA. You would make one set of distribution parameters for each case and run your simulation where each set of distribution parameters is used a fraction $\alpha_i$ of the time.

```math
q(x) = \sum_i \alpha_i q_i(x)
```

Using a mixture of weights is a way to hedge your bets on what about the simulation is causing the outcome to be rare. It does something subtle to the importance calculation because the weight is with respect to the *mixture density.*

```math
w_i = \frac{p(x_i)}{\sum_{k=1}^K \alpha_k q_k(x_i)}
```
In log-space, we would use log-sum-exp.
```math
\log w_i = \log p(x_i) - \log\left(\sum_k\alpha_k q_k(x_i)\right) = \log p(x_i) - \mbox{logsumexp}_k(\log \alpha_k+\log q_k(x_i))
```

Let's say we have three proposal distributions from which we sample evenly, so $\alpha=[1/3, 1/3, 1/3]$. We run each simulation where for each enabling of a clock we pass in a vector of four (4) distributions. The first is the proposal distribution we want to use to generate events for this run. The next three are the actual distribution $p$ and the other proposal distributions.

```julia
using StatsFuns: logsumexp
log_qmix = logsumexp(log.(α) .+ [log_q1, log_q2, log_q3])
log_weight = log_p - log_qmix
w = exp(log_weight)
```

### Adaptive importance sampling (AIS)

Iteratively refit your proposal to minimize weight variance. A typical loop:

 1. Run IS with current proposal.
 2. Reweight samples to estimate sufficient statistics.
 3. Fit a new proposal (e.g. update rate multipliers or shift parameters) using the weighted MLE under the target.

For instance, this could mean adjusting each event's bias factor until the log-weight variance stops decreasing.

### Cross-entropy (CE) method

This structured version of adaptive importance sampling runs
the simulation multiple times. You let it tell you the optimal choice of parameters
to reduce variance in your estimate of the rare event.
One each run, you choose bias parameters $\theta$ to minimize the Kullback-Liebler.

```math
\mbox{KL}(p^*||q_\theta)
```

Here $p^*$ is the conditional distribution on the rare event.

It helps to use common random numbers. Using mixtures is just fine with this
method, and you can select among the mixture proposals. As in other cases, the
self-normalized estimator is more robust than the unbiased estimator.

 1. P.-T. de Boer, D. P. Kroese, S. Mannor, and R. Y. Rubinstein. “A Tutorial on the Cross‑Entropy Method.” *Annals of Operations Research* 134 (2005): 19–67.
 2. R. Y. Rubinstein. “The Cross‑Entropy Method for Combinatorial and Continuous Optimization.” *Methodology and Computing in Applied Probability* 1, no. 2 (1999): 127–190.
 3. R. Y. Rubinstein and D. P. Kroese. *The Cross‑Entropy Method: A Unified Approach to Combinatorial Optimization, Monte‑Carlo Simulation and Machine Learning.* Springer, 2004.

## Variance-reduction companions

### Self-normalized Estimator

Try the self-normalized estimator.

```math
\sum_{i=1}^N \mathcal{I}(x_i>0) w_i / \sum_i w_i
```
The standard estimator.
```math
\sum_{i=1}^N \mathcal{I}(x_i>0) w_i / N
```

### Control variates

If you know some part of the model's behavior exactly (like average promoter activity in our example), then you can leverage it to stabilize the estimate of what you don't know, the rare event probability. The part you know is a *correlated statistic,* $h(x)$, and what you need to know about it is its expectation, $E_p[h]$. Then you can add this into the total expectation without biasing the final estimate.

```math
\hat{\mu}=\sum_i w_i(f(x_i)-c(h(x_i)-E_p[h])).
```

A good choice of control variate has a strong linear correlation with your outcome, such as producing a lot of proteins in our example. It has a known, or easy-to-calculate, expected value. It doesn't have wild swings itself which could amplify noise rather than decrease it.

This technique is often overlooked and can be powerful.

### Stratified or quasi-Monte-Carlo sampling

For path-space simulations, you can stratify on states or numbers of events in the chain, ensuring balanced exploration.
In Julia, you can use QuasiMonteCarlo.jl to drive RNGs for low-discrepancy trajectories.

### Nested or multi-level Importance Sampling

For instance, here we are looking for $P(X>1000)$, so split.

```math
P(X>1000) = P(X>100)P(X>300|X>100)P(X>600|X>300)\cdots
```
This splitting / multilevel approach drastically cuts variance for ultra-rare events.


## Numerical and stability tools

### Log-sum-exp normalization

Use this for any aggregate quantity that involves exp/log across samples.

Do the log-space trick. Just like log-sum-exp in machine learning.
Underflow is a problem. $\exp(-700)$ is small enough to underflow.

```math
\frac{e^{\Delta_i}}{\sum_j e^{\Delta_j}} = \frac{e^{\Delta_i}-\mbox{max}(\Delta)}{\sum_j e^{\Delta_j-\mbox{max}(\Delta)}}
```
This makes the probabilities or expectations identical but improves numerical stability.

```julia
Δ = basal - weighted
Δ = Δ .- maximum(Δ)  # stabilize if you ever aggregate multiple paths at once
importance = exp.(Δ)
```

### Log-domain accumulation

Keep any running-sums of log-likelihoods with log-add-exp rather than exp-sum-log.

### Outlier clipping

In diagnostics (not in the final estimator), cap extermely large weights (say top 0.1%) are recompute ESS. If the estimate hardly changes, your results are stable.

## Interpretability and sanity checks

### Check expected weight

Use log-sum-exp to check that $E_q[w]\approx 1$.

### Check relative contribution

Compute $w_i/\sum_j w_j$ for top samples. If a handful contribute over 50% of the total weight, you're in a weight collapse regime.

### Compare multiple bias directions

Run several small-bias proposals and look for constistent rare-event probability estimates. If results differ wildly, your current proposals are too aggressive.

## For GSMP-specific models

Here each event type has a known distribution, so you can

 * bias only a subset of clocks.
 * conditionally reweight partial paths as soon as one rare component occurs. This can reduce full-path variance.









