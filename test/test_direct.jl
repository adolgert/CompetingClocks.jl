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


@safetestset markov_direct_empty = "MarkovDirect empty hazard" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using ..DirectFixture
    md = MarkovDirect()
    distributions = CTDist[]
    rng = MersenneTwister(90497979)
    current = 0.0
    when, which = next(md, distributions, current, rng)
    @test isinf(when)
    @test which === nothing
end


@safetestset markov_direct_prob = "MarkovDirect probabilities correct" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest, confint
    using ..DirectFixture
    md = MarkovDirect()
    distributions = vcat(
        fill(Exponential(1.5), 10),
        fill(Exponential(1), 10)
    )
    rng = MersenneTwister(90422342)
    hilo = zeros(Int, 2)
    curtime = 0.0
    for i in 1:10000
        when, which = next(md, distributions, curtime, rng)
        hilo[(which - 1) ÷ 10 + 1] += 1
    end
    ci = confint(BinomialTest(hilo[1], sum(hilo), 3 / 5))
    @test ci[1] < 3 / 5 < ci[2]
end
