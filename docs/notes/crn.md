# Common Random Number Design Questions

Common Random Numbers (CRN) are a variance-reduction technique. If your goal is to see how changing the parameters of a simulation changes the mean of the outcome, then first run a simulation while recording the random numbers. Then change parameters but replay random numbers from the first run. This lets you compare means of the outcome with less variance. I haven't seen a lot of mathematical analysis of this technique, and it has its drawbacks. For instance, two simulations may diverge early in a way that makes it impossible to think of reuse the same random draws for the same transitions. Nevertheless, CRN is a common tool that greatly reduces the number of simulation runs required for estimation of the effect of parameter changes.

There are two challenges to implementing CRN. The first is to figure out which draws in a set of simulations might be considered the same in some way. If we are sampling a chemical simulation or a generalized stochastic petri net, we can identify draws by their chemical reactants or by the transition. This isn't foolproof because one simulation may never fire a particular transition while another may fire it many times.

The other challenge is to implement CRN without complicating how the code samples distributions. It would be great to make CRN transparent to the sampling process. It's often implemented as a recording of calls to `rand()`. However, popular distributions, such as Exponential and Erlang, use `randexp()` as their base call, not `rand()`. Julia can make this even more obscure because it has layers of APIs within its sampling system. Let's see what we can work out with two main approaches.


## CRN with Sampling by Inversion

The classic version of CRN is to record, for each clock, for each time it's sampled, the random uniform variate in $[0, 1)$. Then sample the distribution using inversion of the C.D.F. This library has a few samplers that use both sampling of distributions by inversion of the CDF and sampling distributions with the exponential version of the same draw. Here is a snippet that chooses, depending on the distribution, which one to use.

```julia
invert_space(::Type{LinearSampling}, dist, survival) = cquantile(dist, survival)::Float64
invert_space(::Type{LogSampling}, dist, survival) = invlogccdf(dist, survival)::Float64

function sample_by_inversion(
    distribution::UnivariateDistribution, ::Type{S}, te::Float64, when::Float64, survival::Float64
    ) where {S <: SamplingSpaceType}
    if te < when
        te + invert_space(S, truncated(distribution, when - te, Inf), survival)
    else   # te > when
        te + invert_space(S, distribution, survival)
    end
end
```

Because there are two ways to sample distributions, we could implement CRN by storing two different kinds of values, either a uniform draw or an exponential draw. We would need to annotate, clock-by-clock, which one we stored.

```julia
saved_draws = Dict{(Clock, Int), (SamplingSpaceType, Float64)}
```

This would then work for any random number generator.


## CRN with Samplers

The code in Next Reaction doesn't use sampling by inversion when it doesn't have to, because most distributions are more quickly and reliably sampled with other methods. Instead, it uses the `rand()` method for that distribution and random number generator. I expect, in later versions, that we might use Julia's `Distributions.Sampler` interface to use more specific samplers. What you see now is that the code calls `rand()` and then calculates survival, which is generally a forward calculation, not an inversion.

```julia
if te < when
    shifted_distribution = truncated(distribution, when - te, Inf)
    sample = rand(rng, shifted_distribution)
    tau = te + sample
    survival = survival_space(S, shifted_distribution, sample)
else  # te >= when
    # The distribution starts in the future
    sample = rand(rng, distribution)
    tau = te + sample
    survival = survival_space(S, distribution, sample)
end
(tau, survival)
```

Can we support this kind of CRN? Sure, we could save the state of the RNG at each step. This would be nuts for Mersenne Twisters, which have a ton of state, but it's quite reasonable for the Xoshiro RNG, which has four 64-bit values.

```julia
saved_samplers = Dict{(Clock, Int), Sampler}
```

You don't know how many times each clock calls a random number generator when it complets its `rand()`. As a result, it would make sense to skip the RNG forward a bit between calls to `rand()`. It would be possible to instrument an RNG, run a simulation, and check how many times it's called for each draw. Simulations use many more random numbers than you would ever think.

## Design Choices

Let's say we choose a CRN for sampling by inversion. Then we pass run Next Reaction with that CRN. The Next Reaction implementation now has to avoid using `rand()` when it makes samples. It should instead always use inversion to sample.

CRN feels like you're replacing the RNG with something that is no longer random, but it also changes the API to asking for random numbers because each request should also pass in `Clock`, where `Clock` is an identifier. At which point do we pass in the identifier? Right at the moment `rand()` is called on the distribution?

