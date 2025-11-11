```@meta
CurrentModule = CompetingClocks
```

# Sampler Interface

The main interface to the package is the [`SamplingContext`](@ref).
Within that class are the samplers themselves. This describes their interface.

## Constructors

```@docs
CompetingClocks.SamplingContext(::Type{K}, ::Type{T}, rng::R; kwargs...) where {K,T,R<:AbstractRNG}
CompetingClocks.SamplingContext(builder::SamplerBuilder, rng::R) where {R<:AbstractRNG}
```

## Use a Sampler


```@docs
CompetingClocks.enable!(::SSA{K,T}, ::K, ::UnivariateDistribution, ::T, ::T, ::AbstractRNG) where {K,T}
CompetingClocks.reset!(::SSA{K,T}) where {K,T}
CompetingClocks.copy_clocks!(::SSA{K,T}) where {K,T}
CompetingClocks.disable!(::SSA{K,T}, ::K, ::T) where {K,T}
CompetingClocks.next(::SSA{K,T}, ::T, ::AbstractRNG) where {K,T}
```

## Query a Sampler

```@docs
Base.keytype(::SSA{K,T}) where {K,T}
CompetingClocks.timetype(::SSA{K,T}) where {K,T}
CompetingClocks.enabled(::SSA{K,T}) where {K,T}
Base.length(::SSA)
Base.getindex(::SSA{K,T}, ::K) where {K,T}
Base.keys(::SSA)
```
