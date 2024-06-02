```@meta
CurrentModule = CompetingClocks
```

# Interface

These are methods which are defined for any samplers subtyping `<:SSA`, the abstract sampler type.

## Use a Sampler

```@docs
enable!
disable!
next
reset!
copy!
```

## Query a Sampler

```@docs
getindex
keys
keytype
length
```
