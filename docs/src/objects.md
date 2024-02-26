# Objects in Fleck

This describes the kinds of objects in Fleck in order to explain how they work together.


## Distribution

At the heart of simulations are probability distributions. When this code refers to a distribution, it means a [probability distribution function](https://en.wikipedia.org/wiki/Probability_distribution) in time. When we _sample_ a distribution, we randomly draw a time at which an event will happen.

Let's take a Gamma distribution as an example. For the Gamma distribution, the probability distribution function, as you see it in
[Wikipedia's description](https://en.wikipedia.org/wiki/Gamma_distribution),

``f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha - 1}e^{-\beta x}``

For a simulation, we think of the distribution as being over time, and it's over the time since it was first possible for the event to happen, called ``t_e``. We can write the Gamma distribution for ``t-t_e``.

``f(t-t_e) = \frac{\beta^\alpha}{\Gamma(\alpha)}(t-t_e)^{\alpha - 1}e^{-\beta (t-t_e)}``

In Julia, the [univariate distributions](https://juliastats.org/Distributions.jl/stable/univariate/) don't carry their enabling time as a parameter, so we store it separately.


## Competing Clock

Let's say you make a simulation of rabbits eating kibble. A rabbit just picked up a piece of kibble, and the simulation decides there is a Gamma-distributed time at which it will be ready to eat the next bit of kibble. Meanwhile, another rabbit eats the last bit of kibble. That makes it impossible for the first rabbit to eat any kibble, so we say that its event is interrupted. That interruption was a result of the rabbits competing, but what the simulation sees is distributions competing to fire first, and we call them competing clocks.



We have a small naming problem because, while simulations use distributions in time, they sometimes turn them on or off. This combination of a known distribution (such as an Exponential or Gamma distribution) and the ability to temporarily pause it creates its own probability distribution, which we call a clock.

For example, let's make a model of a random walker on a chessboard. Given that the walker is at a grid location, ``(i, j)``, that walker can move to one of four directions, ``(i-1, j)``, ``(i+1, j)``, ``(i, j-1)``, or ``(i, j+1)``. In order to simulate this, we might choose four Exponential distributions, one for each direction. When the walker reaches the side of the chessboard, it doesn't make sense to let it walk off the chessboard, so, at that moment, the simulation _disables_ the ability to walk left. There is an Exponential distribution associated with walking left, but its [hazard rate](https://en.wikipedia.org/wiki/Failure_rate#Failure_rate_in_the_continuous_sense) will be zero while the piece is at the side of the chessboard. The simulation considers "move left" a clock with an Exponential distribution whose firing is disabled.


## Sampler

The main responsibility of Fleck.jl is to provide samplers of competing clocks. Given a list of enabled distributions, decide which one is next to fire and when it will fire.

In order to decide which clock fires next, a sampler needs some information.


### Initialization


The sampler may have a constructor to configure its memory usage or other resources.

### Update Clock State

For a continuous-time simulation, each event happens at a distinct time. Right after an event happens, the simulation is changed, and there are three possible ways clocks are affected.

1. A clock may be disabled. For instance, if the event was crashing a car, you can no longer crash a crashed car.

2. A clock may be enabled. For instance, once a person is infected, they can now infect all neighboring people.

3. The rate of a clock may change. We discussed moves on a chessboard earlier. Maybe there's now a breeze so that the right to move right is greater than the rate to move left. The Exponential distributions would change accordingly.
