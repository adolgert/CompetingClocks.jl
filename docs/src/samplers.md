# Understanding Samplers

From the perspective of someone who uses samplers in a simulation the features are important.

 * Does this sampler work only for Exponential distributions or also non-Exponential distributions?
 * Does it support using atomic distributions, such as delta functions?
 * Does it take advantage of a GPU?
 * Does it put any restrictions on clock keys, such as requirements that they be sequential integers or that there be few keys?
 * Does it support variance reduction techniques, such as
   - Common random numbers
   - antithetic variates
   - importance sampling
   - control variates?
 * Does it maintain accuracy for extreme distribution parameters or extreme draws?

Answers to these questions come from understanding how samplers are built. A modern description of samplers is in Marchetti [Marchetti:2019].

## Series of choices

### Split into marginal and conditional

Fleck supports sampling from generalized semi-Markov processes (GSMP). Every sampler of GSMP is sampling a joint space of which clock event is next, $E$, and which time is next, $T$. The clock event is a discrete choice, and the time is a continuous choice.
The first step to making a GSMP sampler is to split that joint distribution into marginals and conditionals.

The reason to split the joint distribution is that, to sample any joint distribution, the usual method is to sample from the marginal of one random variable and then, given that value, sample the conditional of the other random variable. There are three ways to split a GSMP's joint distribution.

```math
\begin{aligned}
P[E,T] & = P[E]P[T|E] \\
& = P[E|T]P[T] \\
& = min(P[T_i])\quad \forall i
\end{aligned}
```

The first split is rarely, if ever, used because we specify GSMP simulations by the probability in time that each clock will have an event. It takes some work to calculate the marginal over events. If we use notation where $\lambda$ is the hazard, $S$ is the survival, and $t_0$ is the current simulation time, then the equation an integral over a product including all enabled survivals.

```math
P[E] = \int_{t_0}^\infty \lambda_i(t)\prod_i S_i(t) dt
```

Instead, the second split is used all of the time for Exponential distributions. It's the choice that leads to the class of **Direct methods.** That's because the marginal over time is easy to calculate for exponentials.

```math
P[T] = 1 - \exp\left(-\sum_i \int_{t_0}^t \lambda_i(s)ds\right)
```

For exponential distributions, the hazard is constant, which makes the integral trivial. But note that Gillespie's first two papers on exact stochastic simulation included a derivation of this direct method for non-exponential distributions. As with the Exponential distribution, this can be done analytically for the Weibull distribution. For general distributions, the code just has to integrate the right-hand side of the above equation. This is the kind of operation GPUs love.

The third split is a  way to sample the joint distribution by sampling from the conditional distribution of every enabled clock, which we call **competing clocks.** Define a random variate that is the minimum of those clocks, which defines the time to fire. This split matches our instinct that clocks are competing to fire next, and the first one wins.

### Sampling strategy

Now that we've chosen what random variables to sample, in what order, we choose a sampling strategy. By strategy, I mean one of

 * sampling by inversion, also known as [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling)
 * sampling by rejection, also known as an [acceptance-rejection method](https://en.wikipedia.org/wiki/Rejection_sampling)

These are broad methods for sampling, each of which has specific variants for particular distributions or, in the case of GSMP, for random variables defined by sets of competing distributions.

### Specialize for the distribution

Here, we would specialize, for instance, for having all Exponential distributions or all Weibull distributions. Stochastic simulation also introduces another wrinkle in sampling univariate distributions.

In a simulation, it is often the case that a clock is disabled and enabled again, or that a clock that was enabled with one rate is re-enabled with a different rate because it depends on the state of the system, and that state has changed. As a result, we frequently sample clocks that are shifted left. A clock is shifted left when the traditional zero-time for the distribution is in the past. In this case, we have options for how we modify sampling the distribution.

We could use a rejection method, where we sample the distribution and reject samples that are times in the past. We could sample by inversion and select random variates that correspond only to future times. We could sample by rejection and choose a known distribution whose hazard always exceeds that of the shifted distribution. There are lots of options, but none are built into the Julia Distributions package as ways to sample shifted distributions. It's a small complication that applies to stochastic simulation, and we handle it in the code.

### Data structures and algorithms

The next concern is how to express these sampling methods in code. Many of the well-known exact stochastic simulation (SSA) algorithms for chemical simulation are variations, not on how to sample the distribution, but on what data structures and algorithms to use for accelerating that sampling method.

For instance, if we are sampling a simulation using only Exponential distributions, then the Direct method has many named variants. The first step for the Direct method is to sample a time. The time comes from inverting the survival of the system. Here, we denote a random variate between [0,1] by $U$.

```math
U = \exp\left(-\sum_i \lambda_i t\right)
```

The simulation is more effecient if we can maintain $\sum_i \lambda_i$ as a partial sum, so that we can index into it with $\log U$. Further, every time there is a next event in the simulation, the list of enabled hazards changes, so we want to modify that list after each event, while maintaining the sum. For this purpose, there have been many data structures proposed. The general problem is known, in computer science as the [prefix sum or prefix scan](https://en.wikipedia.org/wiki/Prefix_sum). An early attempt was an interesting structure called a [Fenwick tree](https://en.wikipedia.org/wiki/Fenwick_tree). The state of the art for this class of approaches was the [optimized direct method](https://pubmed.ncbi.nlm.nih.gov/15332951/), which is really just one of the prefix scan algorithms applied to the Direct method.

## Examples

### First to fire

The simplest way to sample a GSMP is to sample each clock the moment it's enabled. I don't know that there is a common name for this, so I'm calling it "first to fire." This is a competing clocks split of the joint distribution. For each clock, it can ask the Julia library to sample that clock's distribution in whatever way it sees fit. This is helpful for both speed and accuracy because distributions have methods that are specialized for that distribution and, often, for particular ranges of distribution parameters. For example, it seems easy to sample an exponential using inversion. Here, $U$ is a uniform variate between [0,1], and $t$ is the sample time.

```math
t = (1/\lambda) \log U
```

However, every mathematical library uses the [Ziggurat method](https://en.wikipedia.org/wiki/Ziggurat_algorithm) because it is at least twice as fast. There might be a downside if you want to use automatic differentiation on the code.

### Next reaction method

The next reaction method also samples each clock, but it samples in such a way that it reduces usage of random number generation. Random number generation used to be slow, but it is now no longer a concern. (There's a funny paper about this that I'm having trouble finding.) Nevertheless, there are some reasons to use this method.

When a clock is first enabled, the next reaction method samples a uniform variate in [0,1], and then it finds a putative next event time for that clock using inversion. The original random variate is considered the **total survival for this clock.** The interesting move is that, if that same clock is disabled, this method saves some measure of how far the clock has gotten. That is, it measures the survival of the clock and subtracts that survival from the total survival. If the clock is ever enabled again, its new firing time is determined by the remaining survival.

Is that allowed? The authors of the Next Reaction sort of prove it in ["Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels"](https://pubs.acs.org/doi/full/10.1021/jp993732q), but they can fortunately can rely on the work of Kurtz [Kurtz:1970], which I don't quite see in their references. Nevertheless, Anderson and Kurtz amended this work with ["Continuous time markov chain models for chemical reaction methods"](https://link.springer.com/chapter/10.1007/978-1-4419-6766-4_1).

In Anderson and Kurtz, they take the same approach as the next reaction method. It's still split by the conditional firing times. It's still sampling by inversion, but they change to a log-space for the sampling of individual distributions. That's the whole change. Instead of storing the total survival, $U$, they store the log of that quantity, which turns out to be much more efficient for exponentials and Weibulls, which are used most frequently.

In Fleck, you'll see a single sampler that uses a combination of the next reaction method (Gibson and Bruck) and the modified next reaction method (Anderson and Kurtz). Based on the particular univariate distribution, the code uses a lookup into a performance table to choose either a linear space or log space. It's the best of both worlds, and it's just a matter of changing data structures a little bit.

Why do people like the next reaction method when the first to fire is much less fussy and will use the appropriate sampler for the distribution, every time? It's about features. The next reaction method stores data that helps calculate importance samples. It also makes it easy to implement common random numbers. Finally, the next reaction method always samples from $U=[0,1]$ or from an exponential distribution, and then it transforms that value into the sample of the particular clock's distribution. This is called pathwise sampling, and it enables a simple method for taking derivatives of distributions.

### Direct method for exponentials

The Direct method is really quite solved, as described above. However, there is one complication that has to do with clock keys. Literature about prefix scans assumes that there are integer indices into an ordered, fixed list of integers. That is, you will use the integers 1-100 for the whole simulation. It seems so useful to have clock keys that can be strings, tuples, or other immutable types, so Fleck uses a prefix scan that maintains a dictionary of keys.

If the prefix scan has a dictionary of keys, then it could remember the keys forever or it could forget them. If you have a long-running simulation that uses an ever-growing number of event keys, then it's important to remove keys that are no longer in use. If you have a simulation that uses a fixed set of keys, it's easier to keep them in the dictionary.

It would be interesting to ask whether we could create a data structure that is a keyed prefix sum, instead of using a dictionary that indexes into a prefix sum.

## References

[Kurtz:1970] Kurtz, Thomas G. "Solutions of ordinary differential equations as limits of pure jump Markov processes." Journal of applied Probability 7.1 (1970): 49-58.

[Marchetti:2019] Marchetti, Luca, Corrado Priami, and Vo Hong Thanh. Simulation algorithms for computational systems biology. Vol. 1. Berlin, Germany, 2019.
