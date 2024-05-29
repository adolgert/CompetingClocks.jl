# Next Reaction Method

The reference is ["Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels"](https://pubs.acs.org/doi/10.1021/jp993732q).

In the paper they say that they will speed up the First Reaction Method by doing some things:

  1. store $\tau_{i}$ (putative firing times), not just $a_{i}$ (intensity).
  2. use a dependency graph to only update those two values when needed.
  3. re-use $\tau_{i}$ when appropriate (based on independence of firing times of each clock's Poisson process).
  4. switch from relative time between reactions to absolute time (meaning the $\tau_{i}$ will not have to change if the intensity didn't change).
  5. use efficient data structures to store those two values for each clock.

Call the (directed) dependency graph $G$, where a directed edge from $i$ to $j$ means that the firing of transition $i$ will change the intensity of transition $j$. Then the algorithm is given as:

  1. Initialize:
      - set initial state $X$, set $t\leftarrow 0$, generate $G$.
      - calculate all $a_{i}$.
      - draw putative firing times for each transition distributed as $\tau_{i} \sim \text{Exp}(a_{i})$.
      - store $\tau_{i}$ values in an indexed priority queue $P$.
  2. Let $\mu$ be the transition whose $\tau_{i}$ achieves the minimum.
  3. Set $\tau \leftarrow \tau_{\mu}$
  4. Update state according to $\mu$, and set $t \leftarrow \tau$.
  5. For each transition $\alpha$ that is the target of an edge originating from $\mu$:
      - update $a_{\alpha}$
      - if $\alpha \neq \mu$, set $\tau_{\alpha} \leftarrow (a_{\alpha,\text{old}}/a_{\alpha,\text{new}})(\tau_{\alpha}-t)+t$
      - if $\alpha = \mu$, draw $\rho \sim \text{Exp}(a_{\mu})$ and set $\tau_{\alpha} \leftarrow \rho + t$
      - replace $\tau_{\alpha}$ in $P$ with the new value.
  6. Go to step 2 or quit.

In the general case, for non-Markovian or inhomogeneous Markovian processes, 


During the $n^{\text{th}}$ iteration of the First Reaction method, the random variables giving the relative putative times to firing are $R_{\alpha}$ for each transition. For homogeneous Markov processes we assume they follow an exponential distribution with parameter $a_{\alpha}$. The random variable giving the putative absolute time to firing is $T_{\alpha}=R_{\alpha}+t_{n}$. The distribution of $T_{\alpha}$ is the original shifted to the right by $t_{n}$.

## Implementation in CompetingClocks

Let $\xi_{j}(v)$ be the most recent enabling time of clock $j$. When it is first enabled, we sample by inversion, which solves the following equation for the putative firing time $t^\prime$.

$$
S^{\prime}_{j} = S_{j0}(\xi_{j}(v),t^\prime)
$$

Where $S^{\prime}_{j} = 1 - U_{j}$ and $U_{j} \sim Unif(0,1)$.

If the functional form of the survival changes before the clock fires at some time $t_{1}$ (e.g. due to another clock firing), then we must solve for a new putative firing time (also labeled $t^\prime$). This is called shifting the distribution and it looks like:

$$
S^{\prime}_{j} = S_{j0}(\xi_{j}(v),t_{1}) S_{j1}(t_{1},t^\prime)
$$

If we let the original enabling time be labeled at $t_{0} = \xi_{j}(v)$, and there are $n-1$ changes of functional form that occur (and the clock still has not fired), the general form of shifting the distribution is:

$$
S^{\prime}_{j} = S_{jn}(t_{n},t^\prime) \Pi_{m=0}^{n-1} S_{jm}(t_{m},t_{m+1})
$$

Again, the goal is to solve for $t^\prime$, the updated putative firing time, after the $n-1$-th change.

The MNRM tries to do this efficiently by maintaining the value:

$$
S^{\prime}_{jn} = S^{\prime}_{j} / \Pi_{m=0}^{n-1} S_{jm}(t_{m},t_{m+1})
$$

as computational state to solve:

$$
S^{\prime}_{jn} = S_{jn}(t_{n},t^\prime)
$$

after each change to update $t^\prime$.

## Anderson's method

The modified next reaction method is described here.

Instead of working with survivals, the MNRM works with cumulative hazards.

The general form of "shifting the distribution" for MNRM looks like:

$$
\log(S^{\prime}_{j}) = L_{jn}(t_{n},t^{\prime}_{j}) + \sum_{m=0}^{n-1} L_{jm}(t_{m},t_{m+1})
$$

The MNRM maintains the value of:

$$
L^{\prime}_{jn} = \log(S^{\prime}_{j}) - \sum_{m=0}^{n-1} L_{jm}(t_{m},t_{m+1})
$$

as computational state to solve:

$$
L^{\prime}_{jn} = L_{jn}(t_{n},t^{\prime}_{j})
$$

for the new putative firing time $t^{\prime}_{j}$.