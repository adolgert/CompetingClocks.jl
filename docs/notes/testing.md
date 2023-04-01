# Testing Samplers

## Introduction

We need lots of reassurance that the samplers are correct. From the software side, we can follow the SOCK principle for testability: simple, observable, controllable, knowable. Keep components simple. Give tests a way to enforce invariants. Give tests a way to manipulate components for ease of testing. Ensure that the principles of each component are clearly stated so that they can be tested.

From the statistical side, we have some tools for verifying that samples match known distributions. These are the Kolmogorov-Smirnov test, the Anderson-Darling test, and the multinomial test. These rely on our ability to observe and store empirical distributions.

We can start from two main categories of tests, those that exercise the client interface (and nothing else), and those that look inside the component in some way that requires another method call.


## Risk

We will test most where we are most concerned about invisible failure. These are the failures that don't cause exceptions but do cause incorrect results. Complexity in these samplers comes from the following features:

 1. Distributions can be both non-Markovian.
 2. Distributions can have a starting time that is in the past, immediate, or in the future.
 3. Distributions can be disabled and enabled with different patterns of disabling and enabling.

There are some simple failures for which we should test but that carry less risk. Examples are:

 1. A sampler may work for a clock identifier that is an integer but not for an identifier that is a tuple of a string and an integer.
 2. A sampler may not work for a random number generator that isn't the one used in testing.

The kinds of failures that are more worrisome might include.

 1. A user specifies a non-Exponential distribution for a sampler that only accepts Exponential distributions, but there is no error.
 1. A user specifies an enabling time for a distribution for a sampler that doesn't respect enabling times, because it expects Exponential distributions.
 1. A sampler returns an incorrect sample for distributions that have enabling times in the future because this kind of distribution is more rarely used.
 1. A sampler has an exception when two distributions return the same time, and this error isn't seen in testing because it's a rare occurrence.
 1. A sampler doesn't account properly for a particular pattern of disabling and enabling a transition, and it isn't caught because we account for statistics and fail to check outliers.

It's a worrisome list, but that's why we test!


## Special Cases

Special cases are stochastic models where we know exactly what the result should be. The more complex the system and the more specific the knowledge about the expected result, the more informative the test will be.

Our best example is the [two-state random walk by Weiss](https://github.com/afidd/Semi-Markov/tree/master/doc/manual) (citation below). It's a non-Markovian model with a complicated rule set and an exact marginal.

We know a lot about exact solutions to susceptible-infectious-recovered problems, even some using non-Exponential distributions. We can use some of those.


## Marginals

**Graph Occupancy** - Let's construct an undirected graph and create a transition from each node along the edges. Then run it for some time and measure how long the system spends on each node (the waiting time) and the probability it goes to neighboring nodes. We can run this with a very trustworthy sampler, the First Reaction sampler, and compare with other samples using a multinomial test.

**Walk a Path** - Take a set of distributions and draw a few samples of transitions and times. Find some relatively common path that's a few transitions long. Reject all other paths. Then look at the distribution of transitions out of that state. Do this for several different samplers to see that it's the same. We can vary the types of distributions and any delay of their enabling time. This is both a multinomial test over the probability of next states and a test of the distribution of times.


## Beyond Interface Tests

When we test the external interface, we have know exactly what the input distributions are, and we can observe samples of those distributions. We might test more powerfully by obeying the observability part of the SOCK principle.

We could ask each sampler, at each time step, when every one of the enabled transitions would fire if it were to fire. Even the Next Reaction method knows this value.


## Bibliography

G. H. Weiss, "The two-state random walk," Journal of Statistical Physics, vol 15, no. 2, pp 157-165, 1976.
