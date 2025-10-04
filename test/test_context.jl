# For this macro, need to import it into main for it to work.
import InteractiveUtils: @code_typed

@testset "SamplerContext doesn't penalize missing pieces" begin
    using CompetingClocks
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
    # Assert that when compiled, the if-then statements in the context
    # are compiled away!
    res = @code_typed enable!(sampler, 1, Exponential(0.5), 0.0, 0.0)
    ci = first(res)
    branch_count = count(expr -> isa(expr, Core.GotoIfNot), ci.code)
    @test branch_count == 0
end
