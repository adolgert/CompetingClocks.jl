
# Low-level Sampler Interface

The main interface to the package is the [`SamplingContext`](@ref).
Within that class are the samplers themselves. This describes their interface.

## Common Interface to Low-level Samplers

### Use a Sampler

Low-level sampler interface. Most users should use [`SamplingContext`](@ref) instead,
which is documented in [Context Interface](@ref).

```@docs
CompetingClocks.SSA{Key,Time}
CompetingClocks.enable!(::SSA{K,T}, ::K, ::UnivariateDistribution, ::T, ::T, ::AbstractRNG) where {K,T}
CompetingClocks.reset!(::SSA{K,T}) where {K,T}
CompetingClocks.disable!(::SSA{K,T}, ::K, ::T) where {K,T}
CompetingClocks.next(::SSA{K,T}, ::T, ::AbstractRNG) where {K,T}
```

### Query a Sampler

```@docs
CompetingClocks.timetype(::SSA{K,T}) where {K,T}
CompetingClocks.enabled(::SSA{K,T}) where {K,T}
Base.length(::SSA)
Base.getindex(::SSA{K,T}, ::K) where {K,T}
Base.keys(::SSA)
```

### Copies of Samplers

```@docs
CompetingClocks.clone(sampler::SSA{K,T}) where {K,T}
CompetingClocks.copy_clocks!(::SSA{K,T}) where {K,T}
CompetingClocks.jitter!(sampler::SSA{K,T}, when::T, rng::AbstractRNG) where {K,T}
```

## Individual Samplers

```@docs
CompetingClocks.CombinedNextReaction{K,T}
CompetingClocks.DirectCall{K,T,P}
CompetingClocks.FirstReaction{K,T}
CompetingClocks.FirstToFire{K,T}
CompetingClocks.PSSACR{K,T}
CompetingClocks.Petri{K,T}
CompetingClocks.RSSA{K,T}
```

## Hierarchical Sampling
```@docs
CompetingClocks.MultipleDirect{SamplerKey,K,Time,Chooser}
CompetingClocks.MultiSampler{SamplerKey,Key,Time,Chooser}
```

## Helper Classes

```@docs
CompetingClocks.DebugWatcher{K,T}
CompetingClocks.TrackWatcher{K,T}
CompetingClocks.TrajectoryWatcher{K,T}
CompetingClocks.PathLikelihoods{K,T}
```
