# # Generalized Semi-Markov processes
#
# ## Definition
#
# We said that a semi-Markov process is a set of random variables representing states and times of a system, $\{X_i,T_i\}$. For each state at time $T_i$, the probability of the next state and time is some distribution, $P[X_{i+1},T_{i+1}|X_i,T_i]$. A generalized semi-Markov process (GSMP) also has the same states and times an the same distribution of next states and times, but it's more specific about how to calculate the probability.
#
# For a GSMP, every change in state of the system, from $X_i$ to $X_{i+1}$, is the result of an event $E_j$. Each event is associated with a distribution of event times, also called firing times, and those times are distributed as $f_j(\tau)$ where $\tau=t-T_i$. The probability of the next state and time is determined by the minimum firing time of all events enabled at $(X_i,T_i)$. The distributions of these events are our competing clocks. When one event fires, it changes the state of the system, and, as a result, some events may be disabled and new events may be enabled.
#
# ## Considerations
#
# Once we create this separate object, the event, it raises questions about corner cases. If two events change the state in the same way, then they are the same as one event whose hazard rate is the sum of the two events. If an event does not change the state, then this is like a self-loop in a Markov chain, and it complicates how we count states. As with the semi-Markov process, there remains the question of events which are immediate. If an event can happen at $T_{i+1}=T_i$, then it is possible for a system to fail to progress to a later time, and that's a problem.
#
# The GSMP is a specific form of a semi-Markov process that requires $P[T_{i+1}|X_i,T_i]$ (note the $X_{i+1}$ isn't on the left side) be determined by the minimum time to the next event, so it is determined solely by the $f_j(\tau)$. Each event is also defined by how it changes the state. If events change the state in a deterministic way, such that there is some function $X_{i+1}=\chi_j(X_i)$, then the likelihood of the system is the product of $P[X_{i+1},T_{i+1}|X_i,T_i]$ for each time step, determined by the mimimum of event times. It is traditional to define generalized stochastic Petri nets [Haas:2002] this way, and Anderson and Kurtz's excellent short manuscript presents GSMP in the same light using counting processes [Anderson:2015]. However, Haas describes GSMP as allowing events to be stochastic. This means the likelihood has two terms, $P[E_j,T_{i+1}|X_i,T_i]P[X_{i+1}|E_j, T_{i+1}, X_i,T_i]$.
#
# ## Events and Physical State
#
# Glynn presented GSMP by distinguishing the physical state of the system from the clock state [Glynn:1989], and Shedler is known for having the clearest presentation [Shedler:1987]. He represented the physical state as a set of states $p=(p_1,p_2,p_3...)$. Then each event, $E_j$, is defined in relation to those physical states:
#
#  * An event is *enabled* by an enabling function which depends on a subset of physical states. This is a function $e(\{p_l\},T_i)\rightarrow \text{bool}$. In this notation, the curly bracked $\{\}$ indicate a "set of" something.
#  * When the event is not enabled, it is *disabled.*
#  * The distribution of event times for an event is determined by a subset of physical states. This is a function $f(\{p_m\},T_i)\rightarrow \text{pdf}$.
#  * An event creates a new state by changing some subset of the physical state, and that function can depend on another subset of the physical state, which isn't changed. This is a function $\chi(\{p_j\},\{p_k\},T_i,T_{i+1})\rightarrow \{p'_j\}$.
#
# Look at all the subsets. There is a subset for enabling, transition rates, modified state, and catalyst state (which affects the action but isn't modified). If we think of the physical state as nodes in a graph and the events as nodes in a graph, then each subset associated with an event forms a different kind of edge in a bipartite graph.
#
# The state of the system at any time is more than the physical state. It's the physical state plus the history of when each event clock was enabled. It is even reasonable to include in the state of the system every past event that fired and the time it fired, which is called the *filtration* of the stochastic process.
#
# The idea behind introducing the notion of physical state and subsets of the physical state is to help think about a semi-Markov process where events live longer than a single time step. Glynn wanted to attach those long-lived competing processes to some state because, in practice, there is something about the state of the world that remains the same at each time step. The brilliance of Anderson and Kurtz's monograph is that they start their model as a set of counting processes [Anderson:2015]. Any state of the system is a predictable function of the filtration (event history) of the counting processes.
#
# In the nomenclature of the GSMP, an event defines a change to substates, and every possible pair of states $(X_{i+1},X_i)$ defines a transition. In general, the number of possible transitions is combinatorially larger than the number of possible events, as Haas covers in detail [Haas:2002].
#
# ## Formalisms of GSMP
#
# There are many frameworks for simulation where the simulation is in continuous-time and the next event is determined by competition among clocks.
#
#  * Simulations of chemical reactions. Here, the physical state is chemical components and simulations are usually, but not always, Exponentially-distributed.
#  * Queueing theory models of networks, production, and computation. Here, the state is in queues and the events are reprsented by servers.
#  * Epidemiological models of disease spread among individuals. These models are often hand-coded, but they look a lot like chemical simulation with non-Exponential distributions.
#  * Vector-addition systems are an older form of simulation where the state is a vector of integers and every event gives or takes from the vector of integers. Again, it looks a lot like chemical simulation.
#  * Generalized stochastic Petri nets are what happens when engineers use GSMP. There is a strong vocabulary used to define the state as marking and places, but they conform to what is described above.
#
# ## Extensions to GSMP
#
# ### Atomic Hazards and Differential Equations
#
# Most presentations of GSMP assume that the distribution times of events are continuous distributions that are well-behaved, but the structure of the stochastic process remains well-defined if we allow distributions that have jumps. Examples of such distributions are delta functions. We could for instance, say that an event has a 1/3 chance of happening in 2 minutes an a 2/3 chance of happening in 5 minutes.
#
# More commonly, what if the next step in a simulation were determined by an ordinary differential equation (ODE) that depends on the current state and time? A simulation could, at the enabling time, integrate the ODE to find when it predicts the next event and then enable a distribution that is a delta function centered at that predicted time.
#
# The caveat for atomic hazards it's possible for two atomic hazards to happen at exactly the same time. Fleck doesn't have a way to guarantee which of those events happen first. It certainly doesn't have a way to randomly select which event should happen from a configurable probability distribution. This is a feature specific to samplers that allow instantaneous and simultaneous events. The way we handle the possiblility of simultaneous events is to schedule atomic events relative to continuous-time events or to ensure that there is only one atomic hazard in a whole simulation.
#
# ### Markov Decision Process
#
# In reinforcement learning, there is a model for how the world changes and a model for how to make decisions that depend on the world's history. For a GSMP-based decision process, there are two terms in the likelihood for the next state and time.
#
# ```math
# P[X_{i+1},T_{i+1}|X_i,T_i]P[A_j]
# ```
#
# The additional stochastic variable $A_j$ is the decision at each step.
#
# ### Piecewise-deterministic Markov Process
#
# 
#
# ## References
#
# [Anderson:2015] Anderson, David F., and Thomas G. Kurtz. Stochastic analysis of biochemical systems. Vol. 674. Berlin, Germany: Springer International Publishing, 2015.
#
# [Haas:2002] P. J. Haas, Stochastic Petri Nets: Modelling, Stability, Simulation, Springer-Verlag, New York, New York, USA, 2002.
#
# [Glynn:1989] P. Glynn, A GSMP formalism for discrete event systems, Proceedings of the IEEE 77 (1) (1989) 14– 23, ISSN 00189219, [URL](http://ieeexplore.ieee.org/lpdocs/epic03/wrapper.htm?arnumber=21067).
#
# [Shedler:1987] Shedler, Gerald S. Regeneration and networks of queues. Vol. 3. Springer Science & Business Media, 1987.
#