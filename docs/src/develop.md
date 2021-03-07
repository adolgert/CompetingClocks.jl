# Design

The scope of this is that I want to make a continuous-time simulation that is rooted in stochastic processes and as fast as they can be. The state will be the simplest, a vector of integers. Transitions, however, can be non-Exponential distributions. I'd like to exercise all of the top sampling methods. This code asks how to set up a best-practices continuous-time simulation in Julia. The fastest method for sampling is to use a hierarchical setup, where different sets of transitions get the samplers that are most efficient for their distributions.

The main research question is how to hook the state and transitions of the stochastic process to the sampler algorithm. The mathematical work I've seen does sometimes discuss hierarchical sampling, but it doesn't discuss algorithms for how best to connect a transition with its appropriate sampler. How do we specify it that connection? There's nothing in the math about distributions having a name or an indexing key.

While we normally simulate, there are some other ways to use simulation code, and those other ways give a clue about how to separate functions.

- Fire transitions in a given order without sampling them. Work through a trace.
- Construct a state by specifying an initial state and sequence of transitions.
- Measure likelihood of a state without sampling. This is used for MCMC.

## Vector addition systems

A vector addition system (VAS) defines the state as a vector of non-negative integers. Transitions are two matrices of shape (transition x state). One indicates how large the state vector has to be for that transition to be enabled. Together, they define values removed and added to the state vector when a transition fires.

In practice, vector addition systems aren't implemented with vectors. The professional versions are defined on bipartite graphs where the state is one node type and transitions are the other node type.


## Semi-Markov Vector Addition System

We'll make one little change to the VAS. Instead of having only exponential distributions, it can have general distribution types. That adds complexity to the sampler, but the samplers will handle that. For the VAS, it adds some state and adds burden to specification. The state now needs to track the time at which each transition was enabled, known as ``t_e``. We'll call the main part of the state its physical state, which is the vector of integers in this case.

In order to specify a Markov VAS, you need to give a rate for each transition. A chemical simulation has a very particular function to determine the rate of transitions. For chemicals, the rate is a function of the number of possible combinations of the chemicals. For a Semi-Markov VAS, which we're making, we'll let the state be any function of the state of the system.

There's one other wrinkle to specifying a Semi-Markov VAS. The Zimmerman simulation book points out that non-Exponential transitions can have memory or be memoryless. He means that, if a transition is enabled and disabled, we need to decide whether that affects its rate the next time it's enabled. If the transition were emptying a bathtub, and we put in the stopper plug, then taking out that plug depends on the previous enabling. It must have memory. If it were the time a robot needs to put a gas cap on an assmebly-line car, it would need to begin from scratch during a restart.


## Core responsibilities

1. The state of the system is dictated by stochastic processes theory:
   a. Physical state, which can be a vector of integers for this version.
   b. Enabling times. This is the time at which each transition was, or wasn't enabled.
   c. The last time any transition fired.

2. Sampling depends only on distributions, not any other information about each transition. All optimized samplers keep a cache of the distributions of enabled transitions and modify that cache with each firing.


## Sampling

Most samplers optimize their work by tracking what transitions were last enabled. They treat newly-enabled and newly-disabled transitions as a modification to the likelihood of the next sample. Some samplers don't store any state. First-reaction method is this way, as is the original Direct method. These are naive samplers, and I don't think it's worth optimizing for them. I would, instead, assume that all samplers track which transitions were last enabled. I'll add that as a layer on the first-reaction and direct methods.

## Sequence

The simulation inner loop is critical for efficient simulation. For continuous-time simulation, it's a game of clock-cycles, so we should get this right. This inner loop can have a lot of different kinds of calculations to do, and those calculations benefit from good cache use, right near the processor on L1. As a consequence, we want to design the parts of the inner loop so that they interleave work that operates on the same data.

If we think of a single transition as a consecutive set of steps, each of which has a different function, then those functions traverse state updates, transitions calculations, and sampling decisions.

1. Fire a transition.
2. Modify state.
3. Calculate changes to transitions.
4. Update sampler information.
5. Sample for next transition.

However, the state and transitions will be stored on a graph for the most efficient simulations. We'd like to make the most efficient simulations possible. When state and transitions are stored on a graph, it helps to perform all operations on graph data while doing a _single traversal of the graph._

```
Fire a transition
For each [state changed]
    Update a single state
    For each [transition depending on that state]
        Calculate the hazard rate for that transition.
        Update the sampler about that updated hazard.
```

The goal, then, is to make code where firing, finding affected transitions, and sampling those transitions, is all separate, but they get called nicely in order within a single loop.

## Unwrapping the sequence

I just said the sequence needs to intertwine two algorithms, how we update state and how we sample. Let's start with what the algorithms need to do, so that we get it correct. The complexity comes from doing non-exponential distributions. Each distribution now has an enabling time which moves the time = 0 for that distribution.

The function to update the state and transitions is a moment where we break invariants of the system. That's what makes it difficult. If we change one part of the state and haven't yet changed the transitions that depend on that state, then we are mid-way through what should be an atomic change.

Let's look at the situation without thinking about any caching anywhere. The state of the system is the physical state and the transition state.

* Physical state is a vector of integers, for this Vector Addition System.
* Transition state is a vector of transitions. Each transition has the following.
  1. A univariate distribution over time.
  2. An enabling time to indicate its left bound.
  3. Whether it is enabled (could be encoded in the enabling time, but why.)

Initialization is more interesting than we give it credit. We choose an initial physical state, but we can also choose a transition state that, for instance, includes transitions that were enabled in the past. They would need to be consistent with the current state. We usually start with a boring initialization and run a while in order to do this setup.

Let's send in an event and see what happens. Watch for the event that fires being treated differently.

1. Fire the event. The event is a change in state associated with a transition.
   a. Change time to the event's time.
   b. Mark event's transition as disabled.
   c. Modify the physical state.
   d. Recalculate for every possible transition:
      * Enabled
      * Disabled
      * Changed, meaning its distribution or enabling time is now different.

2. Update sampler with firing, then enabling, disabling, and changing. The same transition that fired can be enabled.

Why did we track the firing transition separately? If we treated it with all other transitions, then the transition might be enabled before the event and still enabled after the event.

## Common representation of transitions

Is there a way that we can create code to help define and use transitions for different simulations?

A transition has these functions:

1. Is it enabled and, if so, with what distribution and time?
2. Fire to change physical state.
3. Fire to change transition state.
