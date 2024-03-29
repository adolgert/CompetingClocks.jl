# Notation for Distributions

This describes implementation of distributions for sampling of
distributions for stochastic simulation in continuous time. This kind of
simulation makes specific demands on what calls a distribution must
support, and those calls are different from what libraries provide. This
is a guide for implementation of new distributions and a way to ensure
that those implemented look correct. If something is wrong here, it
matters, so file a bug report.

## Notation

First, let's affix notation. The cumulative distribution function of
every regular distribution can be written as an integral over its hazard
rate, $\lambda$

$$F(t)=1-e^{-\int_{0}^t \lambda(s)ds}.$$

All algorithms for stochastic simulation treat distributions as being
defined in absolute time, specified as an enabling time, $t_e$,

$$F(t, t_e)=1-e^{-\int_{0}^{t-t_e} \lambda(s)ds}.$$

Working with distributions in absolute time is a simple shift of the
time scale and will be ignored in further discussions, although the
enabling time, $t_e$, will certainly appear in code.

The density function is the derivative of the cumulative distribution
function,

$$f(t)=\frac{dF(t)}{dt}=\lambda(t)e^{-\int_{0}^t \lambda(s)ds}.$$

The survival is

$$G(t)=1-F(t)=e^{-\int_{0}^t \lambda(s)ds}.$$

Because survival is multiplicative, we further label the survival from
time $t_0$ to $t_1$ as

$$G(t_0, t_1)=\frac{G(t_1)}{G(t_0)}=e^{-\int_{t_0}^{t_1} \lambda(s)ds}$$

## Using Julia's Distributions

Julia's continuous univariate distributions support a fixed interface.
In this section, we look at how to translate any distribution into the
operations above.

In this table, `d` is the distribution.

|  Julia call          | Notation                           |
|  ------------------- | ---------------------------------- |
|  `cdf(d,t)`          | $F(t)$                             |
|  `quantile(d,q)`     | $F^{-1}(q)$                        |
|  `logcdf(d,t)`       | $\ln(F(t))$                        |
|  `ccdf(d,t)`         | $G(t)$                             |
|  `logccdf(d,t)`      | $-\int_0^t \lambda(s)ds$           |
|  `quantile(d,q)`     | $F^{-1}(q)$                        |
|  `cquantile(d,q)`    | $F^{-1}(1-q)=G^{-1}(q)$            |
|  `invlogcdf(d,lp)`   | $F^{-1}(e^{l_p})$                  |
|  `invlogccdf(d,lp)`  | $G^{-1}(e^{l_p})$ or $-\int_0^{t(l_p)}\lambda(s)ds=l_p$ |
|  `randexp(rng)`      | $-\ln(1-U)$                        |


## Requirements for a Continuous-Time Simulation

### Shifted Sample

The First Reaction method requires that we sample a distribution given
that we known it has not yet fired by a time $t_0$. The statement that
it hasn't fired by time $t_0$ creates a new distribution from which to
sample. If the old distribution had the hazard $G(t)=G(0, t)$, it could
be written as

$$G(0, t)=G(0, t_0)G(t_0, t).$$

It is this partial survival, since $t_0$, that we want to sample.
Solving for $G(t_0, t)$ and subtracting both sides from 1,

$$1-G(t_0, t)=\frac{G(0, t_0)-G(0, t)}{G(0, t_0)}.$$

Written in terms of the cumulative distribution functions, the cdf of
the new distribution, which we'll call $F(t, t_0)$, is

$$F(t, t_0)=\frac{F(t)-F(t_0)}{1-F(t_0)}$$

This kind of distribution could be sampled by a rejection method, but
the default way to sample it is by inversion, which means generating a
uniform random value between $[0,1]$ and solving $U=F(t)$ for $t$. For
Eq.   `gibson-shifted`{.interpreted-text role="eq"}, this becomes

```math
\begin{aligned}
 U&=F(t,t_0,t_e) \\
  &=\frac{F(t,t_e)-F(t_0,t_e)}{1-F(t_0,t_e)} \\
U(1-F(t_0,t_e))&=F(t,t_e)-F(t_0,t_e) \\
F(t,t_e)&=U(1-F(t_0,t_e))+F(t_0,t_e) \\
F(t-t_e)&=U(1-F(t_0-t_e))+F(t_0-t_e) \\
t-t_e &= F^{-1}\left[U(1-F(t_0-t_e))+F(t_0-t_e)\right] \\
t &= t_e+F^{-1}\left[U(1-F(t_0-t_e))+F(t_0-t_e)\right]
\end{aligned}
```

We will call this operation **SampleShifted.**

### Hazard Rate for Next Reaction

The Next Reaction method requires sampling a distribution such that the
quantile is saved, so that later adjustments to the distribution can use
the same quantile.

During a simulation, the hazard rate, $\lambda$, is a function of the
state of the system, $X(t)$. The state of the system only changes in
jumps, so the hazard rate is effectively a series of distributions in
time. For instance, a hazard rate, from enabling time $T_0$ to firing
time $T_3$, might have three parts.

```math
\begin{aligned}
  \lambda(\{X_0, T_0\}, t)=h_0(t) & \qquad T_0 \le t < T_1 \\
  \lambda(\{X_0, T_0, X_1, T_1\}, t)=h_1(t) & \qquad T_1 \le t < T_2 \\
  \lambda(\{X_0, T_0, X_1, T_1, X_2, T_2\}, t)=h_2(t) & \qquad T_2 \le t < T_3
\end{aligned}
```

The algorithm therefore samples for a firing time from $h_0(t)$ as soon
as the transition is enabled, but that time will turn out to be wrong
(we call it a putative time). Later, the algorithm will resample using
$h_1(t)$ using the original sample's quantile and taking time off the
clock. If the first sample were by inversion, it would look like solving
this equation for $t$ (still ignoring enabling times),

$$U=1-\exp\left(\int_0^{t}h_0(s)ds\right).$$

Then a later sample would use the same $U$, but with knowledge that the
distribution now contains a new part, $h_1(t)$,

```math
\begin{equation}\tag{cdf-consume}
U=1-\exp\left(-\int_0^{t_1}h_0(s)ds\right)\exp\left(-\int_{t_1}^{t}h_1(s)ds\right)
\end{equation}
```

Anderson had the bright idea to write the quantile as an exponential
quantile, $1-U=e^{\ln (1-U)}$, so that the equation requires only
addition of integrated hazards,

```math
\begin{aligned}
  \ln(1-U)&=-\int_0^{t_1}h_0(s)ds-\int_{t_1}^{t}h_1(s)ds \\
  \int_{t_1}^{t}h_1(s)ds&=-\ln(1-U)-\int_0^{t_1}h_0(s)ds.
\end{aligned}
```

As the underlying distribution, $h_i(t)$, changes, the right hand side
gets smaller and smaller. Let's call the sum of all consumed hazard
$\gamma$,

$$\gamma=\sum_i \int_{t_i}^{t_{i+1}}h_i(s)ds$$

The algorithm therefore needs three operations from the distribution.

1.  **MeasuredSample**\-\--Sample the distribution, returning the
    exponential quantile. Calling the random number generator, `rng`
    and the putative time $t_p$, it's

    $$\left(rng, h_0(t)\right) \mapsto \left(t_p, -\ln(1-U)\right).$$

2.  **ConsumeSample**\-\--Consume remaining quantile for the next
    sample. If the sum of hazard in the past is called $\gamma$, then

    $$\left(\gamma, h_i(t), t_i\right) \mapsto \left(\gamma', t_{i+1}\right)$$

3.  **Putative**\-\--Generate a new putative time from the exponential
    quantile and the consumed hazard,

    $$\left(-\ln(1-U), \gamma, t_i, h_i(t)\right) \mapsto p_t$$

The nice part about the first step is that there is no need to sample by
inversion. Any sampling method will do, as long as the exponential
quantile is calculated.

### Cumulative Distributions for Next Reaction

The original form of the Next Reaction, by Gibson and Bruck, was written
in terms, not of the hazards, but of the cumulative distribution
functions. This form remains useful because some distributions are much
simpler, or more accurate, to sample as cdfs instead of sampling from
their hazard rates.

Returning to Eq. cdf-consume above, this can be rewritten as

$$1-U=\exp\left(-\int_0^{t_1}h_0(s)ds\right)\exp\left(-\int_{t_1}^{t}h_1(s)ds\right)$$

In terms of the survival functions, this becomes

$$1-U=G_0(0, t_1)G_1(t_1, t)$$

If we wish to solve this for $t$, then, in terms of the survival, it
looks like

$$G_1(t_1, t)=\frac{1-U}{G_0(0, t_1)}$$

Writing the left-hand side as a cumulative distribution function
requires the transformation

$$G(a, b)=\frac{G(b)}{G(a)}=\frac{1-F(b)}{G(a)}$$

so we have

$$F_1(t)=1-G_1(t_1) \frac{1-U}{G_0(0, t_1)}$$

This generalizes to many steps as

$$F_j(t)=1-G_j(t_j) (1-U) \prod_i^j \frac{G_i(t_i)}{G_i(t_{i+1})}$$

Let's call the running product on the right $\delta$,

$$\delta=\prod_i^j \frac{G_i(t_i)}{G_i(t_{i+1})}$$

Then the algorithm requires three operations

1.  **MeasuredSample**\-\--Sample the distribution, returning the
    quantile. Calling the random number generator, "rng," and the
    putative time $t_p$, it's

    $$\left(\rng, h_0(t)\right) \mapsto \left(t_p, 1-U\right).$$

2.  **ConsumeSample**\-\--Consume remaining quantile for the next
    sample.

    $$\left(\delta_i, h_i(t), t_i\right) \mapsto \left(\delta_{i+1}, t_{i+1}\right)$$

3.  **Putative**\-\--Generate a new putative time from the exponential
    quantile and the consumed hazard,

    $$\left(1-U, \delta_i, t_i, h_i(t)\right) \mapsto p_t$$

As you can see by comparison with the hazards version, it's simple to
write the algorithm to accommodate either method of sampling. Therefore,
each distribution can choose which interface to support.
