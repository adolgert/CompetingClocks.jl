```@meta
CurrentModule = CompetingClocks
```

# Interface

These are methods which are defined for any samplers subtyping `<:SSA`, the abstract sampler type.

```@docs
SamplingContext
```

## Use a Sampler

```@docs
enable!
disable!
next
reset!
fire!
copy_clocks!
```

## Query a Sampler

```@docs
getindex
keys
keytype
timetype
length
```
