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