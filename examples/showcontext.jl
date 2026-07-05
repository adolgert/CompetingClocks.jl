using CompetingClocks
using InteractiveUtils
using Random
using Distributions
K = Int64
T = Float64
S = FirstToFire{Int64,Float64}
# Type params: {K, T, Sampler, RNG, Watchers<:Tuple, CRN, DelayedState}
# No watchers, no CRN, no delayed state => empty watcher tuple and Nothing.
SC = CompetingClocks.SamplingContext{K,T,S,Xoshiro,Tuple{},Nothing,Nothing}

rng = Xoshiro(899987987)
# Fields: sampler, rng, watchers, crn, split_weight, time, fixed_start,
#         sample_distribution, delayed
sampler = SC(FirstToFire{Int64,Float64}(), rng, (), nothing, 1.0, 0.0, 0.0, 1, nothing)
for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
    enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0)
end
ci = first(@code_typed enable!(sampler, 1, Exponential(0.5), 0.0, 0.0))

branch_count = count(expr -> isa(expr, Core.GotoIfNot), ci.code)
@assert branch_count == 0
