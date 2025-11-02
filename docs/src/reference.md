```@meta
CurrentModule = CompetingClocks
```

# Samplers

The choice of sampler determines specific algorithms that are used to sample, update, and disable clocks. Helpers also exist that are useful for logging, utilizing common random numbers, and hierarchical sampling.

## Sampler Supertype

```@docs
SSA
SamplerBuilder
add_group!
```

## Sampler Types

```@docs
FirstReaction
FirstToFire
DirectCall
CombinedNextReaction
```

## Sampling Helpers

```@docs
CommonRandom
freeze_crn!
reset_crn!
split!
haskey
misscount
misses
enabled
MultiSampler
ChatReaction
Petri
DebugWatcher
TrackWatcher
consume_survival
sampling_space
steploglikelihood
FromInclusion
```
