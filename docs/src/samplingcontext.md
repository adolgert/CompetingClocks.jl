# Sampling with SamplingContext

This is the main interface to the package.

 1. **Build the SamplingContext.** The [`SamplerBuilder`](@ref) helps choose what hierarchical
    samplers to use within the [`SamplingContext`](@ref). You will see the builder
    gives you options for debugging output, calculation of likelihoods,
    use of variance reduction (common random numbers), and more.
    

```julia
ClockKey = Int64
Time = Float64
rng = Xoshiro(899987987)
builder = SamplerBuilder(ClockKey, Time)
sampling = SamplingContext(builder, rng)
```

You can configure the builder by adding specific samplers for different kinds
of clock keys. Each set of keys specified in an [`add_group!`](@ref) goes to
its own sampler in the hierarchy.

 2. **Write a sampling loop**. The sampler will tell you what event comes
    [`next`](@ref) and then you choose that event using [`fire!`](@ref).
    However, you write your own state for the system and update that state
    according to the event. Firing the event disables that event.

```julia
function one_epoch(model, sampler)
    observe = Observer()
    step_model!(model, sampler, (:on, 0), time(sampler))
    when, which = next(sampler)
    while !isnothing(which)
        fire!(sampler, which, when)
        step_model!(model, sampler, which, when)
        when, which = next(sampler)
        observer(model, when, which)
    end
    return observer
end
```

 3. **Enable/disable events during state update.** When updating the model
    state, you might enable or disable events. Here is an example from
    gene expression.

```julia
function translate_protein(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.mrna)
    model.protein += 1
    if pre_event_total > 0
        transrate1 = Exponential(inv(θ[:translate][1] * pre_event_total))
        transrate2 = Exponential(inv(θ[:translate][2] * pre_event_total))
        enable!(sampler, (:translate, 0), [transrate1, transrate2])
    end
end
```
