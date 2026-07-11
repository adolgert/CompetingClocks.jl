# Debugging a Simulation that Uses CompetingClocks

It can be difficult to debug simulations that use time-varying hazards, especially
when it's possible to cancel events and re-enable events with new hazard rates.
This library provides a few tools to help.

## Sample Unlikely Events with the Petri Sampler

If your simulation has rare events, and you want to test code paths for those
rare events, one approach is to sample the next clock entirely randomly. The
Petri sampler ignores the distribution of times for events and picks the
next event evenly among all enabled keys.

```julia
using CompetingClocks
using Distributions
using Random

K = Int
T = Float64
petri_sampler = CompetingClocks.Petri{K,T}(seed=90210)

now = 0.0
enable!(petri_sampler, 3, Exponential(100.0), now, now)
when, what = next(petri_sampler, now)
@assert what == 3
@assert abs(when - now - 1) <= 1e-9
```

The Petri sampler is a low-level sampler, so it is reached through a qualified
name, `CompetingClocks.Petri`, and called with the low-level `(sampler, clock,
distribution, enabling_time, when)` form of `enable!`. The sampler owns its own
random streams, selected by the seed, so no random number generator is passed.
Because it ignores each clock's distribution, the time it returns is simply one
plus the current time.


## Keep History of Enabling and Disabling

The `DebugWatcher` records every time an event is enabled or disabled.
You wouldn't normally want to use all of this time and memory (and memory churn),
but it can be helpful to see the trace of all events enabled and disabled,
in addition to those that fired.

```julia
for enabling_event in enabled_history(sampler)
    println("$(enabling_event.clock), $(enabling_event.when)")
end
for disabling_event in disabled_history(sampler)
    println("$(disabling_event.clock), $(disabling_event.when)")
end
```

Look at [`CompetingClocks.enabled_history`](@ref) and [`CompetingClocks.disabled_history`](@ref).
