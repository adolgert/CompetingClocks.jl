using SafeTestsets


@safetestset first_reaction_smoke = "FirstReaction smoke" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: MersenneTwister
    using Distributions: Exponential, Gamma

    rng = MersenneTwister(90422342)
    min_when = Inf
    seen = Set{Int}()
    sample_time = 0.5
    for i in 1:100
        sampler = FirstReaction{Int}()
        enable!(sampler, 1, Exponential(1.7), 0.0, 0.0, rng)
        enable!(sampler, 2, Gamma(9, 0.5), 0.0, 0.0, rng)
        enable!(sampler, 3, Gamma(2, 2.0), 0.0, 0.0, rng)
        disable!(sampler, 2, sample_time)
        when, which = next(sampler, sample_time, rng)
        push!(seen, which)
        @test when > sample_time
        min_when = min(min_when, when)
    end
    @test min_when < sample_time + 0.01
    @test seen == Set([1, 3])
end


@safetestset first_reaction_empty = "FirstReaction empty" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: MersenneTwister

    rng = MersenneTwister(90422342)
    sampler = FirstReaction{Int}()
    when, which = next(sampler, 5.7, rng)
    @test when == Inf
    @test which === nothing
end


@safetestset first_reaction_binomial = "FirstReaction binomial" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: MersenneTwister
    using ..DirectFixture: test_exponential_binomial

    rng = MersenneTwister(12349678)
    sampler = FirstReaction{Int}()
    test_exponential_binomial(sampler, rng)
end


@safetestset first_reaction_weibull_binomial = "FirstReaction weibull binomial" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: MersenneTwister
    using ..DirectFixture: test_weibull_binomial

    rng = MersenneTwister(12967847)
    sampler = FirstReaction{Int}()
    test_weibull_binomial(sampler, rng)
end


@safetestset first_reaction_single = "FirstReaction single transition" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: Xoshiro
    using HypothesisTests: pvalue, ExactOneSampleKSTest
    using Distributions: Gamma, Exponential, Weibull, truncated
    rng = Xoshiro(8367109004)
    rand(rng, 100)  # burn some early numbers

    sampler = FirstReaction{Int}()
    dist = Weibull()
    sample_cnt = 1000
    enable!(sampler, 1, dist, 0.0, 0.0, rng)
    samples = [next(sampler, 0.0, rng)[1] for i in 1:sample_cnt]
    @test all(isfinite(x) for x in samples)
    @test all(x > 0 for x in samples)
    ks_test = ExactOneSampleKSTest(samples, dist)
    @test pvalue(ks_test) > 0.04
end


@safetestset first_reaction_single_later = "FirstReaction single transition later" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: Xoshiro
    using HypothesisTests: pvalue, ExactOneSampleKSTest, confint
    using Distributions: Gamma, Exponential, Weibull, truncated
    rng = Xoshiro(8367109004)
    rand(rng, 100)  # burn some early numbers

    sampler = FirstReaction{Int}()
    dist = Weibull()
    sample_cnt = 1000
    enable!(sampler, 1, dist, 0.0, 0.0, rng)
    later = 0.7
    samples = [next(sampler, later, rng)[1] for i in 1:sample_cnt]
    @test all(isfinite(x) for x in samples)
    @test all(x > 0 for x in samples)
    ks1_test = ExactOneSampleKSTest(samples, dist)
    @test pvalue(ks1_test) < 0.04
    ks2_test = ExactOneSampleKSTest(samples, truncated(dist, later, Inf))
    @test pvalue(ks2_test; tail = :both) > 0.04
end


@safetestset first_reaction_single_future = "FirstReaction single future transition" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: Xoshiro
    using HypothesisTests: pvalue, ExactOneSampleKSTest, confint
    using Distributions: Gamma, Exponential, Weibull, truncated
    rng = Xoshiro(8367109004)
    rand(rng, 100)  # burn some early numbers

    sampler = FirstReaction{Int}()
    dist = Weibull()
    future = 2.7
    sample_cnt = 1000
    enable!(sampler, 1, dist, future, 0.0, rng)
    later = 0.7
    samples = [next(sampler, later, rng)[1] for i in 1:sample_cnt]
    @test all(isfinite(x) for x in samples)
    @test all(x > future for x in samples)
    ks1_test = ExactOneSampleKSTest(samples, dist)
    @test pvalue(ks1_test) < 0.04
    shifted = [(x - future) for x in samples]
    ks2_test = ExactOneSampleKSTest(shifted, dist)
    @test pvalue(ks2_test; tail = :both) > 0.04
end
