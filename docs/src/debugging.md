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
petri_sampler = Petri{K,T}()
enable!(tw, 3, Exponential(100.0), 0.0, 0.0, rng)
when, what = next(petri_sampler, now, rng)
@assert abs(when - now - 1) <= 1e-9
```

It also returns a time that is one plus the current time.


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
