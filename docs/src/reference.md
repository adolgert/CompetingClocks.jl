```@meta
CurrentModule = CompetingClocks
```

# Samplers

The choice of sampler determines specific algorithms that are used to sample, update, and disable clocks. Helpers also exist that are useful for logging, utilizing common random numbers, and hierarchical sampling.


## Sampling Helpers

```@docs
CommonRandom
freeze_crn!
reset_crn!
split!
misscount
misses
ChatReaction
consume_survival
sampling_space
steploglikelihood
FromInclusion
```
