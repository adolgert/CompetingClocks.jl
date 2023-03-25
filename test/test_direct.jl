using SafeTestsets

module DirectFixture
using Distributions: ContinuousUnivariateDistribution
const CTDist = ContinuousUnivariateDistribution
# Use import so that we can extend the method.
using Random: AbstractRNG


export CTDist
end


@safetestset direct_call_prob = "DirectCall probabilities" begin
    using Fleck: DirectCall, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential

    dc = DirectCall{Int}()
    rng = MersenneTwister(90422342)
    propensities = [0.3, 0.2, 0.7, 0.001, 0.25]
    for (i, p) in enumerate(propensities)
        enable!(dc, i, Exponential(p), 0.0, 0.0, rng)
    end
    when, which = next(dc, 100.0, rng)
    @test when > 100
    @test 1 <= which
    @test which <= length(propensities)
end
