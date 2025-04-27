using SafeTestsets


@safetestset direct_call_prob = "DirectCall probabilities" begin
    using CompetingClocks: DirectCall, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential

    dc = DirectCall{Int,Float64}()
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

@safetestset DirectCall_interface = "DirectCall basic interface" begin
    using CompetingClocks
    using Random: Xoshiro
    using Distributions

    sampler = DirectCall{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[1]
    @test keytype(sampler) <: Int64

    @test_throws ErrorException enable!(sampler, 1, Dirac(1.0), 0.0, 0.0, rng)

    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        enable!(sampler, clock, Exponential(when_fire), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[1] == 1/7.9

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "1")

    disable!(sampler, 1, 0.0)

    @test_throws KeyError sampler[1]
    @test sampler[2] == 1/12.3

end

@safetestset direct_call_empty = "DirectCall empty hazard" begin
    using CompetingClocks: DirectCall, next, enable!, reset!
    using Random: MersenneTwister
    md = DirectCall{Int,Float64}()
    rng = MersenneTwister(90497979)
    current = 0.0
    when, which = next(md, current, rng)
    @test isinf(when)
    @test which === nothing
    reset!(md)
end


@safetestset direct_call_later = "DirectCall probabilities correct at a later time" begin
    using CompetingClocks: DirectCall, next, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest, confint
    rng = MersenneTwister(2343979)

    # Given 10 slow distributions and 10 fast, we can figure out
    # that the marginal probability of a low vs a high is 1 / (1 + 1.5) = 3/5.
    # Check that we get the correct marginal probability.
    md = DirectCall{Int,Float64}()
    for i in 1:10
        enable!(md, i, Exponential(1), 0.0, 0.0, rng)
    end
    for i in 11:20
        enable!(md, i, Exponential(1.5), 0.0, 0.0, rng)
    end
    hilo = zeros(Int, 2)
    curtime = 2.5
    for i in 1:10000
        when, which = next(md, curtime, rng)
        hilo[(which - 1) ÷ 10 + 1] += 1
    end
    ci = confint(BinomialTest(hilo[1], sum(hilo), 3 / 5))
    @test ci[1] < 3 / 5 < ci[2]
end


@safetestset direct_call_prob = "DirectCall probabilities correct" begin
    using CompetingClocks: DirectCall, next, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential
    using HypothesisTests: BinomialTest, confint
    using ..DirectFixture: test_exponential_binomial
    rng = MersenneTwister(223497123)
    md = DirectCall{Int,Float64}()
    test_exponential_binomial(md, rng)
end


@safetestset direct_call_copy = "DirectCall copy" begin
    using CompetingClocks: DirectCall, enable!, next
    using Random: MersenneTwister
    using Distributions: Exponential

    src = DirectCall{Int,Float64}()
    dst = DirectCall{Int,Float64}()
    rng = MersenneTwister(90422342)
    enable!(src, 1, Exponential(), 0.0, 0.0, rng)
    enable!(src, 2, Exponential(), 0.0, 0.0, rng)
    enable!(dst, 3, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy!(dst, src)
    @test length(dst) == 2
    enable!(src, 5, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    enable!(dst, 6, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3
end
