using Test
using Distributions
using Random
using QuadGK
using CompetingClocks: cumulative_hazard


@testset "TranscriptionRate Distribution" begin

    @testset "Constructor" begin
        # Valid construction
        d = TranscriptionRate(2.0, 0.5; t0=1.0)
        @test d.α_max == 2.0
        @test d.k_rem == 0.5
        @test d.t0 == 1.0

        # Default t0
        d_default = TranscriptionRate(2.0, 0.5)
        @test d_default.t0 == 0.0

        # Invalid parameters
        @test_throws ErrorException TranscriptionRate(-1.0, 0.5)
        @test_throws ErrorException TranscriptionRate(2.0, -0.5)
        @test_throws ErrorException TranscriptionRate(2.0, 0.5; t0=-1.0)
    end

    @testset "Parameter extraction" begin
        d = TranscriptionRate(3.0, 0.8; t0=2.5)
        p = params(d)
        @test p == (3.0, 0.8, 2.5)
    end

    @testset "CDF properties" begin
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # CDF is 0 at t < 0
        @test cdf(d, -0.1) == 0.0

        # CDF is 0 at t=0, then positive for t > 0
        @test cdf(d, 0.0) == 0.0
        @test cdf(d, 0.01) > 0.0

        # CDF is monotonically increasing
        @test 0.0 < cdf(d, 0.5) < cdf(d, 2.0) < cdf(d, 5.0) < cdf(d, 10.0) < 1.0

        # CDF approaches 1 as t → ∞
        @test cdf(d, 100.0) > 0.99
    end

    @testset "PDF properties" for t0 in [0.0, 0.5, 1.0, 2.0]
        d = TranscriptionRate(2.0, 0.5; t0=t0)

        # PDF is 0 for t < 0
        @test pdf(d, -0.1) == 0.0

        # PDF is positive for t ≥ 0 (forward shift means hazard starts higher)
        @test pdf(d, 0.0) ≥ 0.0
        @test pdf(d, 0.5) > 0.0
        @test pdf(d, 1.0) > 0.0
        @test pdf(d, 2.0) > 0.0
        @test pdf(d, 5.0) > 0.0
        @test 1.0 - quadgk(x -> pdf(d, x), 0, Inf)[1] < 1e-6
    end

    @testset "logPDF properties" begin
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # logPDF is -∞ for t < 0
        @test logpdf(d, -0.1) == -Inf

        # logPDF is finite for t ≥ 0
        @test isfinite(logpdf(d, 0.0))
        @test isfinite(logpdf(d, 0.5))
        @test isfinite(logpdf(d, 1.0))
        @test isfinite(logpdf(d, 2.0))

        # logPDF matches log(PDF)
        for t in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
            @test logpdf(d, t) ≈ log(pdf(d, t)) atol = 1e-10
        end
    end

    @testset "Numerical stability near t=0" begin
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # Test very close to t=0 (where we start sampling)
        small_eps = [1e-10, 1e-8, 1e-6, 1e-4, 1e-2]
        for eps in small_eps
            t = eps

            # logpdf should be finite (not NaN or -Inf for t ≥ 0)
            lp = logpdf(d, t)
            @test !isnan(lp)
            @test isfinite(lp)

            # pdf should be non-negative and finite
            p = pdf(d, t)
            @test p >= 0.0
            @test isfinite(p)

            # CDF should be non-negative and finite
            c = cdf(d, t)
            @test c >= 0.0
            @test isfinite(c)
        end
    end

    @testset "Cumulative hazard function" begin
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # Cumulative hazard is 0 at t=0
        @test cumulative_hazard(d, 0.0) == 0.0

        # Cumulative hazard is 0 for t < 0
        @test cumulative_hazard(d, -0.5) == 0.0

        # Cumulative hazard is positive for t > 0
        @test cumulative_hazard(d, 0.5) > 0.0
        @test cumulative_hazard(d, 1.0) > 0.0

        # Cumulative hazard is monotonically increasing
        Λ0 = cumulative_hazard(d, 0.5)
        Λ1 = cumulative_hazard(d, 2.0)
        Λ2 = cumulative_hazard(d, 5.0)
        Λ3 = cumulative_hazard(d, 10.0)
        @test 0.0 < Λ0 < Λ1 < Λ2 < Λ3

        # Verify relationship: S(t) = exp(-Λ(t))
        for t in [0.5, 2.0, 5.0, 10.0]
            Λ = cumulative_hazard(d, t)
            S_from_Λ = exp(-Λ)
            S_from_cdf = 1.0 - cdf(d, t)
            @test S_from_Λ ≈ S_from_cdf atol = 1e-10
        end
    end

    @testset "Random sampling" begin
        Random.seed!(42)
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # Generate samples
        n_samples = 1000
        samples = [rand(d) for _ in 1:n_samples]

        # All samples should be >= 0 (support is [0, ∞))
        @test all(s >= 0.0 for s in samples)

        # Samples should have reasonable spread
        @test minimum(samples) >= 0.0
        @test maximum(samples) > 1.0

        # Check that empirical CDF roughly matches theoretical CDF
        test_points = [0.5, 1.0, 2.0]
        for t in test_points
            empirical_cdf = count(s <= t for s in samples) / n_samples
            theoretical_cdf = cdf(d, t)
            # Allow for sampling variability
            @test empirical_cdf ≈ theoretical_cdf atol = 0.1
        end
    end

    @testset "Different parameter regimes" begin
        # Small k_rem (slow approach to α_max)
        d_slow = TranscriptionRate(1.0, 0.1; t0=0.0)
        @test pdf(d_slow, 0.5) < pdf(d_slow, 5.0)

        # Large k_rem (fast approach to α_max)
        d_fast = TranscriptionRate(1.0, 10.0; t0=0.0)
        @test pdf(d_fast, 0.5) > pdf(d_fast, 5.0)

        # Large α_max
        d_high_rate = TranscriptionRate(100.0, 1.0; t0=0.0)
        @test cdf(d_high_rate, 0.5) > cdf(TranscriptionRate(1.0, 1.0; t0=0.0), 0.5)
    end

    @testset "Consistency between PDF, CDF, and quantiles" begin
        d = TranscriptionRate(2.0, 0.5; t0=1.0)

        # d/dt CDF(t) ≈ PDF(t)
        t = 5.0
        h = 1e-6
        numerical_derivative = (cdf(d, t + h) - cdf(d, t - h)) / (2h)
        @test numerical_derivative ≈ pdf(d, t) rtol = 1e-4
    end

    @testset "Edge case: t0 = 0 (unshifted)" begin
        d = TranscriptionRate(2.0, 0.5; t0=0.0)

        # Should behave like the base distribution
        @test cdf(d, 0.0) == 0.0
        @test pdf(d, 0.0) == 0.0
        @test pdf(d, 1.0) > 0.0

        # Random sampling should work
        Random.seed!(123)
        samples = [rand(d) for _ in 1:100]
        @test all(s >= 0.0 for s in samples)
    end

    @testset "Sequential transcription events" begin
        # Simulate the intended use case: sequential transcription events
        Random.seed!(789)
        α_max = 2.0
        k_rem = 0.5

        # First transcription starts at t0=0
        d1 = TranscriptionRate(α_max, k_rem; t0=0.0)
        t1 = rand(d1)
        @test t1 >= 0.0

        # Second transcription starts with t0=t1 (continuing the hazard ramp-up)
        d2 = TranscriptionRate(α_max, k_rem; t0=t1)
        Δt2 = rand(d2)  # Time from t1 to t2
        t2 = t1 + Δt2
        @test Δt2 >= 0.0

        # Third transcription
        d3 = TranscriptionRate(α_max, k_rem; t0=t2)
        Δt3 = rand(d3)
        t3 = t2 + Δt3
        @test Δt3 >= 0.0

        # Verify that the hazard rate increases across events
        # At the start of each waiting period, hazard should be higher
        λ1_start = α_max * (1 - exp(-k_rem * 0.0))  # λ(0) for first event
        λ2_start = α_max * (1 - exp(-k_rem * t1))   # λ(0) for second event = λ(t1) from base
        λ3_start = α_max * (1 - exp(-k_rem * t2))   # λ(0) for third event = λ(t2) from base

        @test λ1_start < λ2_start < λ3_start
        @test λ1_start == 0.0  # First event starts with no transcription activity
        @test λ2_start > 0.0   # Second event starts with some activity
        @test λ3_start > λ2_start  # Third event starts with even more activity
    end
end
