# Context Interface

## Building a Context

```@docs
CompetingClocks.SamplingContext(::Type{K}, ::Type{T}, rng::R; kwargs...) where {K,T,R<:AbstractRNG}
CompetingClocks.SamplingContext(builder::SamplerBuilder, rng::R) where {R<:AbstractRNG}
CompetingClocks.SamplerBuilder
CompetingClocks.add_group!
```

## Choosing Sampler Methods

```@docs
CompetingClocks.NextReactionMethod
CompetingClocks.DirectMethod
CompetingClocks.FirstReactionMethod
CompetingClocks.FirstToFireMethod
CompetingClocks.PartialPropensityMethod
CompetingClocks.PetriMethod
CompetingClocks.RejectionMethod
```


## Using a Context

```@docs
CompetingClocks.enable!
CompetingClocks.disable!
CompetingClocks.next
CompetingClocks.fire!
CompetingClocks.reset!
CompetingClocks.copy_clocks!
Base.time(::SamplingContext)
CompetingClocks.sample_from_distribution!
CompetingClocks.clone
CompetingClocks.enabled
CompetingClocks.enabled_history
CompetingClocks.disabled_history
CompetingClocks.EnablingEntry
CompetingClocks.DisablingEntry
Base.length(::SamplingContext)
CompetingClocks.isenabled
Base.keytype(::SamplingContext)
CompetingClocks.timetype
CompetingClocks.pathloglikelihood
```
