```@meta
CurrentModule = CompetingClocks
```

# Algorithms

Many samplers depend on data structures to allow efficient querying of clocks ordered with respect to some value, usually the firing time. These types and methods implement them for CompetingClocks.

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
