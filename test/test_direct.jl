using SafeTestsets


@safetestset markov_direct_blah = "MarkovDirect initial" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    md = MarkovDirect()
    distributions = fill(Exponential(1.5), 10)
    rng = MersenneTwister(90422342)
    when, which = next(md, distributions, rng)
    @test when > 0
    @test 1 <= which
    @test which <= 10
end


@safetestset markov_direct_empty = "MarkovDirect empty hazard" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    md = MarkovDirect()
    distributions = CTDist[]
    rng = MersenneTwister(90497979)
    when, which = next(md, distributions, rng)
    @test isinf(when)
    @test which === nothing
end


@safetestset markov_direct_prob = "MarkovDirect probabilities correct" begin
    using Fleck: MarkovDirect, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest
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
