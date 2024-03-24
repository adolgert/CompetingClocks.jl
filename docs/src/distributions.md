# How Fleck Decides What Happens Next

Fleck simulates generalized semi-Markov processes (GSMP). Every event in a generalized semi-Markov process is chosen as the result of a competion among clocks to see which fires next. Let's describe here how clocks compete and begin by contrasting that with how SimJulia works.

The flow of control in a [SimJulia](https://simjuliajl.readthedocs.io/en/stable/welcome.html) simulation is based on tasks. Each task performs some action on the state, rolls the dice, and sets a wake-up time. It might wake up as expected, or it might be interrupted by another task's actions. The simulation code decides when it will wake next, and then it draws random numbers to determine what to do. In contrast, a simulation using Fleck will create a set of possible next events, assign a probability distribution for *when* each can happen, and the timing of which happens first determines *which* next event happens. Let's look at how a probability distribution describes the time for an event to happen and then how they compete in Fleck.

## Distributions in Time

Let's say you have a cold. You know you aren't going to recover immediately, but, as days go by, you're more and more sure you'll recover soon. You can read this graph below as the probability to recover given that you have not yet recovered.

![](assets/gammahazard.png)

It starts at zero, meaning there's no way you'll recover as soon as you're sick. It gets more likely over time that you're at the tail end of being sick. This is called a *hazard rate,* and the hazard rate shown is that of a Gamma distribution, commonly used to describe the rate of recovery for a population of individuals who are sick.

If, instead, you want to see the number of people who recover on any given day, that is called a probability distribution function.

![](assets/gammapdf.png)

Where the hazard rate is an instantaneous quantity at a point in time, the probability distribution function (pdf) integrates over all possible future times. If we call the hazard rate $\lambda(t)$ and call the pdf $f$, we get this relationship.

```math
f(t) = \lambda(t) e^{-\int_0^t \lambda(s)ds}
```

The graph of the pdf tells us that the most likely time for this event is a little before time 5, in whatever units. You will see graphs of [pdfs on Wikipedia](https://en.wikipedia.org/wiki/Gamma_distribution#/media/File:Gamma_distribution_pdf.svg) because this is how people usually think about the probability an event happens at some time.

A simulation, however, has multiple events possible at any one time. One event may happen, and then other events need to restart. Let's ask, if you still have a cold on day 5, what is the probability distribution function for when you will recover?

![](assets/shiftedgamma.png)

The probability distribution function changes now that you know you didn't recover earlier than day 5. On the other hand, the hazard rate for recovery from the cold will be unchanged. Using the same hazard rate, we can recalculate the pdf from the new time $t_0=5$.

```math
f(t;t>t_0) = \lambda(t) e^{-\int_{t_0}^t \lambda(s)ds}
```

In some sense, the underlying nature of when an event will happen is described more by the hazard rate.

## Competition

### Individual Distributions

Let's think of a moment when there are three possible next events. There is a Gamma distribution for when you recover from a cold, a Weibull distribution for when you decide to take medicine for the cold, and an Exponential distribution for when your Mom calls you. Each one is described by a distribution in time, and we can think of them as three hazard rates.

![](assets/individual_distributions.png)

The separate hazard rates are what we put into the simulation. Given their competition, the hazard rates will remain unchanged, but the pdfs will change.

### Total Probability

Given three events determined by the three clock distributions above, there is a marginal probability that any one of them will be the first to fire. We calculate this with an integral over the distribution of each event, multiplied by the survivals of the other events.

```math
P[E_i] = \int_0^\infty f_i(t) \prod_{j\ne i} S_j(t) dt
```

That gives the graph on the left.

![](assets/conditional_on_which.png)

The graph on the right shows the conditional distribution in time for each event, given that it was the one that fired, so it is $P[t_i | E_i]$. You can see that these distributions don't match the distributions for the individual events. They are modified by competition.


### Total Time

From the graph above, if we pick a time, $t_1=10$, we can read from the graph three hazard rates, $(\lambda_1(t_1),\lambda_2(t_1), \lambda_3(t_1))$. Each hazard rate is the rate, per unit time, of that event. We know that, if the simulation makes it to $t=10$ without any event happening, the conditional probability for any one of those events is the ratio of hazard rates.

```math
P[E_1|t=t_1] = \frac{\lambda_1(t_1)}{\lambda_1(t_1)+\lambda_2(t_1)+\lambda_3(t_1)}
```

If all distributions are Exponential, then they all have a constant hazard rate, and that conditional probability stays the same over time. If any distributions are non-Exponential, then that conditional probability changes.

![](assets/conditional_on_when.png)

On the left of this graph is the pdf for the first event of the three to fire. We can see this as a marginal $P[t]$ and then the right-hand graph as the conditional $P[E_i|t]$.

## Specification of a Simulation

We simulated using competing clocks when we measure using survival analysis. One is the inverse of the other.

A classic example of survival analysis is a hospital study where participants enter a trial for a drug, and doctors observe the date at which each patient recovers. How would you draw a probability distribution function for such a trial? You would draw a histogram where each week shows a count of how many recovered. What about the hazard rate for such a trial? Here, you would look only at the number of participants remaining in the trial and ask what fraction recover each week. That's a rate of recovery given that participants have not yet recovered.

There are many disciplines that use survival analysis to measure rates of processes. Epidemiology measures rates of contagion, death, and recovery. Demography measures time to have children or move to another country. Reliability studies the time for a part to break and be fixed. Drug discovery measures the time to approve a chemical for the next trial stage. For traffic flow, it's the time to the next intersection. For package delivery, it's the time to the next package scan. All of these do use, or can use, survival analysis to observe a system to simulate.
