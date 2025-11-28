using CompetingClocks
using InteractiveUtils
using Random
using Distributions
K = Int64
T = Float64
S = FirstToFire{Int64,Float64}
SC = CompetingClocks.SamplingContext{K,T,S,Xoshiro,Nothing,Nothing,Nothing}

rng = Xoshiro(899987987)
sampler = SC(FirstToFire{Int64,Float64}(), rng, nothing, nothing, nothing, 0.0)
for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
    enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0)
end
ci = first(@code_typed enable!(sampler, 1, Exponential(0.5), 0.0, 0.0))

branch_count = count(expr -> isa(expr, Core.GotoIfNot), ci.code)
@assert branch_count == 0
