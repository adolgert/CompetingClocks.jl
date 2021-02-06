using Distributions
using HypothesisTests
using Random
using SafeTestsets
using Test


const CTDist = ContinuousUnivariateDistribution

function hazards(callback::Function, distlist::Array{T, 1}, rng) where T <: CTDist
    for dist_idx in 1:length(distlist)
        callback(dist_idx, distlist[dist_idx], 0.0, true, rng)
    end
end


@safetestset markov_direct_blah = "MarkovDirect initial" begin
    md = MarkovDirect()
    distributions = fill(Exponential(1.5), 10)
    rng = MersenneTwister(90422342)
    when, which = next(md, distributions, rng)
    @test when > 0
    @test 1 <= which
    @test which <= 10
end


@safetestset markov_direct_empty = "MarkovDirect empty hazard" begin
    md = MarkovDirect()
    distributions = CTDist[]
    rng = MersenneTwister(90497979)
    when, which = next(md, distributions, rng)
    @test isinf(when)
    @test which === nothing
end


@safetestset markov_direct_prob = "MarkovDirect probabilities correct" begin
    md = MarkovDirect()
    distributions = vcat(
        fill(Exponential(1.5), 10),
        fill(Exponential(1), 10)
    )
    rng = MersenneTwister(90422342)
    hilo = zeros(Int, 2)
    for i in 1:10000
        when, which = next(md, distributions, rng)
        hilo[(which - 1) ÷ 10 + 1] += 1
    end
    ci = confint(BinomialTest(hilo[1], sum(hilo), 3 / 5))
    @test ci[1] < 3 / 5 < ci[2]
end
