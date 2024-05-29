# # Hierarchical Samplers
#
# ## Overview
#
# Continuous-time samplers have to process each event separately. As simulations grow in size, samplers use more memory, and they take more time to select the next event. One approach to speed up sampling is to use multiple samplers arranged in a hierarchy.
#
# ## Why Hierarchical Samplers
#
# There are two reasons to use hierarchical samplers.
#
# Some samplers are better for some simulations. If all transitions are Exponentially-distributed, then an optimized Direct sampler can be the fastest. If all of the distributions are of Exponential families, then Anderson's method is faster than the Next Reaction method. You can split the events among samplers according to the sampler that best fits the behavior of those events.
#
# The other reason to use multiple samplers has to do with the frequency and locality of the events, in the same way we think of the frequency and locality of memory accesses for cache use. If there is a small subset of events that regenerate frequently, it can make sense even to use a [`FirstReaction`](@ref) sampler for those events. While FirstReaction doesn't use a complicated data structure to optimize, it can be winningly fast for small numbers of events. Or, for a spatial simulation, you could make separate samplers for separate parts of the landscape, so that each event tends to affect a limited number of samplers.
#
# ## First Sampler to Fire
#
# All of the different samplers find the first event to fire. If we set up two samplers, so that each holds mutually distinct enabled event distributions, we can ask each sampler which it thinks will fire next. The first to fire is the first of the two samplers. This generalizes to any number of samplers. We can make a [`MultiSampler`](@ref) which contains multiple samplers and always returns the soonest of those contained.
#
# Even further a MultiSampler can contain a MultiSampler if that makes sense.
#
# ## Multiple Direct Samplers
#
# It's possible to make a Direct-style sampler that is hierarchical, too. A Direct sampler works in two steps. It sums the hazards of all enabled events and then selects a time according to the sum of the hazards. The main algorithm of a Direct sampler is to sum hazards. A hierarchical Direct algorithm sums the sums of the hazards and then selects a time.
#
# While hierarchical samplers can contain multiple-Direct samplers, multiple-Direct samplers can only contain other multiple-Direct samplers.
#
# ## How to Split a Simulation
#
# Each time a simulation calls `enable!` and `disable!` for an event, it specifies a key for that event. If the sampler is hierarchical, it can use that key, and maybe the type of the distribution, to choose which sampler handles any given event.
#
# Fleck's approach in the [`MultiSampler`](@ref) type is to let the user specify a function that takes as input the key and the distribution and returns some ID for the chosen sampler. The MultiSampler then remembers that choice for this event key, from that point on.
#
