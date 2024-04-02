```@meta
CurrentModule = Fleck
```

# Reference

 * [Interface](@ref)
 * [Samplers](@ref)
 * [Algorithms](@ref)


## Interface

These are methods which are defined for any samplers subtyping `<:SSA`, the abstract sampler type.

### Use a Sampler

```@docs
enable!
disable!
next
reset!
```

### Query a Sampler

```@docs
getindex
keys
keytype
length
```

## Samplers

The choice of sampler determines specific algorithms that are used to sample, update, and disable clocks. Helpers also exist that are useful for logging, utilizing common random numbers, and hierarchical sampling.

### Sampler Supertype

```@docs
SSA
```

### Sampler Types

```@docs
FirstReaction
FirstToFire
DirectCall
CombinedNextReaction
consume_survival
sampling_space
```

### Sampling Helpers

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
```

## Algorithms

Many samplers depend on data structures to allow efficient querying of clocks ordered with respect to some value, usually the firing time. These types and methods implement them for Fleck.

```@docs
CumSumPrefixSearch
BinaryTreePrefixSearch
KeyedKeepPrefixSearch
KeyedRemovalPrefixSearch
choose
setindex!
rand
set_multiple!
```
