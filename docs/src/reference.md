```@meta
CurrentModule = Fleck
```

# Reference

 * [Interface](@ref)
 * [Samplers](@ref)
 * [Algorithms](@ref)


## Interface


### Use a Sampler

```@docs
enable!
disable!
next
```

### Query a Sampler

```@docs
getindex
keys
keytype
length
```

## Samplers

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
ChatReaction
CommonRandomRecorder
reset!
misscount
misses
DebugWatcher
MultiSampler
SingleSampler
TrackWatcher
```

## Algorithms

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
