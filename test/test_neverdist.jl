using CompetingClocks

@safetestset never_fair = "Never distribution" begin
    using CompetingClocks: Never
    using Distributions: params, partype, mean, median, mode, var, skewness, kurtosis
    using Distributions: pdf, logpdf, cdf, ccdf, quantile, mgf, cf
    using Random

    @testset "construction" begin
        n = Never()
        @test isa(n, Never{Float64})
        @test params(n) == ()
        @test partype(n) == Float64

        n32 = Never{Float32}()
        @test isa(n32, Never{Float32})
        @test partype(n32) == Float32
    end

    @testset "moments and statistics" begin
        n = Never()
        @test mean(n) == typemax(Float64)
        @test median(n) == typemax(Float64)
        @test mode(n) == typemax(Float64)
        @test var(n) == typemax(Float64)
        @test skewness(n) == 0.0
        @test kurtosis(n) == 0.0
    end

    @testset "probability functions" begin
        n = Never()
        # pdf is 0 everywhere - no probability mass at any finite point
        @test pdf(n, 0.0) == 0.0
        @test pdf(n, 1.0) == 0.0
        @test pdf(n, 1e100) == 0.0

        # logpdf is -Inf everywhere
        @test logpdf(n, 0.0) == -Inf
        @test logpdf(n, 1.0) == -Inf

        # cdf is 0 everywhere - event never occurs by any finite time
        @test cdf(n, 0.0) == 0.0
        @test cdf(n, 1e100) == 0.0

        # ccdf is 1 everywhere - survival is certain
        @test ccdf(n, 0.0) == 1.0
        @test ccdf(n, 1e100) == 1.0

        # quantile is always typemax - any quantile maps to "infinity"
        @test quantile(n, 0.0) == typemax(Float64)
        @test quantile(n, 0.5) == typemax(Float64)
        @test quantile(n, 1.0) == typemax(Float64)
    end

    @testset "transforms" begin
        n = Never()
        # mgf: M(t) = E[e^{tX}] for X=∞
        # t <= 0: e^{t·∞} = 0
        # t > 0: e^{t·∞} = ∞
        @test mgf(n, -1.0) == 0.0
        @test mgf(n, 0.0) == 0.0
        @test mgf(n, 1.0) == Inf

        # cf is undefined for X=∞, returning 0 as sentinel
        @test cf(n, 1.0) == 0.0
    end

    @testset "random sampling" begin
        n = Never()
        rng = Random.MersenneTwister(42)

        # Single sample
        @test rand(rng, n) == typemax(Float64)

        # Array sampling
        arr = zeros(5)
        Random.rand!(rng, n, arr)
        @test all(x -> x == typemax(Float64), arr)
    end

    @testset "type stability" begin
        n32 = Never{Float32}()
        @test mean(n32) == typemax(Float32)
        @test typeof(mean(n32)) == Float32
        @test skewness(n32) == 0.0f0
        @test typeof(skewness(n32)) == Float32
    end
end
