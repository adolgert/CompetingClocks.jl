```@meta
CurrentModule = Fleck
```

# Samplers

The choice of sampler determines specific algorithms that are used to sample, update, and disable clocks. Helpers also exist that are useful for logging, utilizing common random numbers, and hierarchical sampling.

## Sampler Supertype

```@docs
SSA
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
CommonRandomRecorder
freeze
misscount
misses
MultiSampler
SingleSampler
ChatReaction
DebugWatcher
TrackWatcher
consume_survival
sampling_space
```
