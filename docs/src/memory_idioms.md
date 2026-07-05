# Re-enabling and Memory

When a clock's situation changes during a simulation — a rate depends on state,
and the state just moved — you have to decide what happens to the time that clock
has *already* spent enabled. For an exponential clock the answer never matters:
the exponential distribution is memoryless, so its future is independent of its
past. For every other distribution it matters a great deal, and it is a
**modeling choice**, not a detail the library can make for you.

The choice enters through the enabling time. At the context layer,

```julia
enable!(ctx, clock, dist, relative_te)
```

sets the clock's zero-time to `te = time(ctx) + relative_te`. The underlying
sampler works in absolute time, where the same call is

```julia
enable!(sampler, clock, distribution, enablingtime, when, rng)
```

with `enablingtime` the absolute zero-time of `distribution` and `when` the
current simulation time. When `enablingtime == when` (equivalently
`relative_te == 0`) the clock starts fresh from *now*. Shifting `enablingtime`
into the past is what encodes memory. There are three idioms.

## Idiom 1: re-enable fresh

Call `enable!` again with a new distribution and `relative_te == 0` (the
default). The clock forgets whatever time it had accumulated and starts over from
the current moment:

```julia
enable!(ctx, clock, Weibull(2.0, 3.0))   # relative_te defaults to 0
```

This is the right model when the change *resets* the process — the physical
mechanism that was counting down has been replaced by a new one. A machine that
is repaired and put back into service, a queue that is emptied and restarted: the
age that had accumulated is no longer meaningful.

## Idiom 2: left-shift via the enabling time

Re-enable with the enabling time in the past, `relative_te < 0`:

```julia
age = time(ctx) - original_enabling_time
enable!(ctx, clock, Weibull(2.0, 3.0), -age)
```

Now the distribution's zero is `age` units before now, so the clock *remembers*
how long it has already been enabled. This samples the distribution conditioned
on the event not having fired before the current time — exactly the shifted,
conditional distribution derived in [Shifted Sample](shifting.md). Use this when
the underlying process keeps running through the state change: the component keeps
aging, the customer keeps waiting, and only the future hazard is being updated.

## Idiom 3: explicit conditioning with `truncated`

You can state the conditioning in the distribution itself, using `truncated`
from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/truncate/),
while keeping the left shift in the enabling time:

```julia
using Distributions
age = time(ctx) - original_enabling_time
enable!(ctx, clock, truncated(Weibull(2.0, 3.0); lower=age), -age)
```

The truncation makes the conditioning explicit — this clock has already survived
to `age` — and the shifted enabling time keeps the distribution's zero at the
true enabling moment, so firing times are measured correctly. This is equivalent
to Idiom 2 and is the clearest choice when the conditioning is more elaborate
than a single left shift — for instance, conditioning on both a lower and an
upper bound.

Beware: truncation *alone*, with `relative_te == 0`, is a different (and usually
wrong) model. The truncated sample is then measured from *now*, so every firing
lands at least `age` in the future, rather than being the remaining life of a
clock that is already `age` old.

!!! warning "These are different stochastic processes"
    Idiom 1 and Idioms 2–3 do not describe the same model. Re-enabling fresh
    throws away the clock's age; left-shifting (or truncating) preserves it. They
    give different firing-time distributions and therefore different sample paths
    and different likelihoods. Because exponential clocks hide the distinction,
    picking the wrong idiom by accident is the most common *silent* modeling error
    in non-exponential simulation. Decide deliberately, per clock, whether a
    change of state should reset a clock or let it keep aging.

The `enable!` docstring spells out the sign convention: with the low-level
signature `enable!(sampler, clock, distribution, enablingtime, when, rng)`,
choosing `enablingtime > when` delays the event, `enablingtime < when` shifts it
left (memory), and `enablingtime == when` starts it fresh. At the context layer
the same three cases are `relative_te > 0`, `relative_te < 0`, and
`relative_te == 0`.

## Fire versus disable

Two verbs end a clock's enabled episode, and they are not synonyms.

Call `fire!` for the clock that *happened*. Firing realizes the draw: the
firing time has been observed, nothing about that draw remains unknown, and if
the clock is enabled again later it starts from a fresh draw.

Call `disable!` for clocks whose *preconditions vanished* — the event that
fired removed their reason to run. Disabling censors the draw: all you learned
is that the clock had not yet fired, so a sampler is entitled to keep the
draw's remaining randomness and reuse it when the clock is enabled again.
`CombinedNextReaction` does exactly that; it is why the next-reaction method
needs only one random number per clock per enabling episode, and it is what
keeps common random numbers aligned across runs.

For exponential clocks the distinction is invisible — memorylessness makes a
fresh draw and a reused draw the same in law and in cost. For everything else,
using the wrong verb is a statistics bug, not a style choice:

!!! warning "Use the verb that matches what happened"
    Calling `disable!` on a clock that actually fired tells a reusing sampler
    to keep residual randomness for a draw that was fully consumed — a later
    re-enable can resurrect a dead draw. Calling `fire!` on a clock that was
    merely switched off throws away residual randomness it was entitled to
    keep — the law of your process survives, but draw reuse and common random
    numbers do not.

The high-level interface makes this hard to get wrong: `fire!(ctx, clock,
when)` on the [`SamplingContext`](@ref) does the right thing for the fired
clock, and your model's `disable!` calls should only ever be about
preconditions. If you drive a sampler directly through the low-level API, the
choice of verb is yours on every event.

Like the older question of whether an event always changes state or only
sometimes does, memory is part of what *you* mean by an event. The library keeps
the clocks consistent; it is up to you to say which process you are simulating.
