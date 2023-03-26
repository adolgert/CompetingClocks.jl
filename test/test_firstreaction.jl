using SafeTestsets

module DirectFixture
using Distributions: ContinuousUnivariateDistribution
const CTDist = ContinuousUnivariateDistribution
# Use import so that we can extend the method.
import Fleck: hazards
using Random: AbstractRNG

function hazards(callback::Function, distlist::Array{T, 1}, rng::AbstractRNG) where T <: CTDist
    for dist_idx in 1:length(distlist)
        callback(dist_idx, distlist[dist_idx], true)
    end
end
export CTDist
end


@safetestset markov_direct_initial = "MarkovDirect initial" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using ..DirectFixture
    md = MarkovDirect()
    distributions = fill(Exponential(1.5), 10)
    rng = MersenneTwister(90422342)
    current = 0.0
    when, which = next(md, distributions, current, rng)
    @test when > 0
    @test 1 <= which
    @test which <= 10
end
