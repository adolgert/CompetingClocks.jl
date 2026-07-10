```@meta
CurrentModule = CompetingClocks
```

# Samplers

The choice of sampler determines specific algorithms that are used to sample, update, and disable clocks. Helpers also exist that are useful for logging, common random numbers via keyed streams, and hierarchical sampling.


## Sampling Helpers

Common random numbers are provided by each sampler's per-clock keyed streams
rather than a recorder: two samplers built from the same seed draw identically
per clock. See [`KeyedStreams`](@ref), [`rekey_streams!`](@ref), and
[`clone`](@ref).

```@docs
KeyedStreams
rekey_streams!
similar_sampler
split!
consume_survival
sampling_space
steploglikelihood
FromInclusion
```
