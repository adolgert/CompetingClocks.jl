# Vector Addition System

One way to understand how to use CompetingClocks is to look at a very simple simulation. A [vector addition system](https://en.wikipedia.org/wiki/Vector_addition_system) (VAS) is a lot like chemical simulations, but it's free of some of the assumptions about chemical rates. The physical state of the system is a vector of integers. Such a simple physical state can make it easier to understand possible complications in how we define events.

There are variations on how to define a vector addition system. Let's begin with one and introduce variations as they become interesting.

A vector addition system is

 * a physical state $p$ which is a set of $d$ integers, labeled $(p_1, p_2,\cdots, p_d)$.

 * a system time, $t$

 * a finite set of events $E$ where each event has an enabling rule, a distribution, and a firing rule.

   - a rule for when the event is enabled, where the rule is an invariant on the physical state, such as $r_i(p)\ge 0$. When the invariant is met, the event is enabled. When the invariant is not met, the event is disabled. For a vector addition system, this rule must be expressible as a matrix multiplication, so $M\cdot \vec{p} .\ge 0$, where the $.\ge$ means we are checking elementwise that every element is non-negative.

   - a distribution in time, where the distribution is determined at enabling time and is a function of the physical state and system time at enabling. The distribution in time can have shocks, such as a delta function, but it may not have a shock at time zero, must have a measure-zero probability of firing at the same moment it is enabled.

   - a rule for how the event changes physical state when it happens. The domain and co-domain of the rule are physical states of the system. Let's limit this model to require all events to modify the state, so that we exclude events that don't modify the state. As with enabling, this must be expressible as a matrix operation on the physical states to produce a new physical state.

   From the rules, we see that an event's state is either disabled or enabled at an enabling time $t_e$.

How would we implement this system? It's set up so that enabling rules and firing rules can be done with linear algebra. We can write this in Julia code 



Simulation using the VAS requires concrete implementation of the interface:

```
zero_state(vas::VectorAdditionSystem)
vas_delta(vas::VectorAdditionSystem, transition_idx)
vas_initial(vas::VectorAdditionSystem, initial_state)
fire!(visitor, vas::VectorAdditionSystem, state, modify_state, rng)
simstep!(fsm::VectorAdditionFSM, state_update::Function, rng::AbstractRNG)
```

The method `fire!` first modifies state, then, for each clock in the system, it checks whether it has been newly enabled or disabled.

  - newly disabled: `visitor` is called to disable that transition
  - newly enabled: `visitor` is called to enable that transition and cache its newly calculated intensity
  - still enabled: we check if the new intensity differs from the old and if so use `visitor` to update the intensity

The method `simstep!` first applies `fire!`, followed by `next`.
