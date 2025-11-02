# Context Interface

## Building a Context

```@docs
SamplerBuilder
add_group!
```

## Using a Context

```@docs
SamplingContext
CompetingClocks.enable!(::SamplingContext{K,T}, ::K, ::Any, ::T) where {K,T}
CompetingClocks.enable!(::SamplingContext{K,T}, ::K, ::Any) where {K,T}
CompetingClocks.enable!(::SamplingContext{K,T}, ::K, ::Vector, ::T) where {K,T}
CompetingClocks.enable!(::SamplingContext{K,T}, ::K, ::Vector) where {K,T}
CompetingClocks.disable!(::SamplingContext{K,T}, ::K) where {K,T}
CompetingClocks.fire!(::SamplingContext{K,T}, ::K, ::T) where {K,T}
CompetingClocks.reset!(::SamplingContext{K,T}) where {K,T}
CompetingClocks.copy_clocks!(::SamplingContext{K,T}, ::SamplingContext{K,T}) where {K,T}
CompetingClocks.split!(::AbstractVector{S}, ::SamplingContext{K,T}) where {K,T,S<:SamplingContext}
CompetingClocks.time(::SamplingContext)
CompetingClocks.sample_from_distribution!(::SamplingContext, ::Any)
CompetingClocks.clone(::SamplingContext{K,T,Sampler,RNG,Like,CRN,Dbg}) where {K,T,Sampler,RNG,Like,CRN,Dbg}
CompetingClocks.freeze_crn!(::SamplingContext)
CompetingClocks.reset_crn!(::SamplingContext)
CompetingClocks.enabled(::SamplingContext)
Base.length(::SamplingContext)
CompetingClocks.isenabled(::SamplingContext{K}, ::K) where {K}
Base.keytype(::SamplingContext{K}) where {K}
CompetingClocks.timetype(::SamplingContext{K,T}) where {K,T}
CompetingClocks.steploglikelihood(::SamplingContext, ::Any, ::Any)
CompetingClocks.pathloglikelihood(::SamplingContext, ::Any)
```
