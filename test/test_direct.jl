using SafeTestsets


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


@safetestset direct_call_empty = "DirectCall empty hazard" begin
    using Fleck: DirectCall, next, enable!
    using Random: MersenneTwister
    md = DirectCall{Int}()
    rng = MersenneTwister(90497979)
    current = 0.0
    when, which = next(md, current, rng)
    @test isinf(when)
    @test which === nothing
end


@safetestset direct_call_prob = "DirectCall probabilities correct" begin
    using Fleck: DirectCall, next, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest, confint
    rng = MersenneTwister(2343979)

    # Given 10 slow distributions and 10 fast, we can figure out
    # that the marginal probability of a low vs a high is 1 / (1 + 1.5) = 3/5.
    # Check that we get the correct marginal probability.
    md = DirectCall{Int}()
    for i in 1:10
        enable!(md, i, Exponential(1), 0.0, 0.0, rng)
    end
    for i in 11:20
        enable!(md, i, Exponential(1.5), 0.0, 0.0, rng)
    end
    hilo = zeros(Int, 2)
    curtime = 0.0
    for i in 1:10000
        when, which = next(md, curtime, rng)
        hilo[(which - 1) ÷ 10 + 1] += 1
    end
    ci = confint(BinomialTest(hilo[1], sum(hilo), 3 / 5))
    @test ci[1] < 3 / 5 < ci[2]
end


@safetestset direct_call_later = "DirectCall probabilities correct at later time" begin
    using Fleck: DirectCall, next, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest, confint
    rng = MersenneTwister(223497123)

    # Given 10 slow distributions and 10 fast, we can figure out
    # that the marginal probability of a low vs a high is 1 / (1 + 1.5) = 3/5.
    # Check that we get the correct marginal probability.
    md = DirectCall{Int}()
    for i in 1:10
        enable!(md, i, Exponential(1), 0.0, 0.0, rng)
    end
    for i in 11:20
        enable!(md, i, Exponential(1.5), 0.0, 0.0, rng)
    end
    hilo = zeros(Int, 2)
    curtime = 20.3
    for i in 1:10000
        when, which = next(md, curtime, rng)
        hilo[(which - 1) ÷ 10 + 1] += 1
    end
    ci = confint(BinomialTest(hilo[1], sum(hilo), 3 / 5))
    @test ci[1] < 3 / 5 < ci[2]
end
