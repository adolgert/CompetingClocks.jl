using Test
using CompetingClocks: hazard
using Distributions
using ForwardDiff

@testset "hazard function" begin
    @testset "Exponential distribution" begin
        d = Exponential(2.0)  # rate = 1/2 = 0.5
        te = 1.0

        # Before enabling time: hazard is 0
        @test hazard(d, te, 0.5) == 0.0
        @test hazard(d, te, te - eps()) == 0.0

        # At and after enabling time: hazard = rate
        @test hazard(d, te, te) == rate(d)
        @test hazard(d, te, 5.0) == rate(d)
    end

    @testset "general UnivariateDistribution (Weibull)" begin
        # Weibull with shape k and scale s has hazard h(t) = (k/s) * (t/s)^(k-1)
        k = 2.0  # shape
        scale = 3.0
        d = Weibull(k, scale)
        te = 1.0

        # Before enabling time: hazard is 0
        @test hazard(d, te, 0.5) == 0.0
        @test hazard(d, te, te - eps()) == 0.0

        # After enabling time: hazard = pdf / ccdf
        for t in [te, 2.0, 5.0, 10.0]
            tau = t - te
            expected = pdf(d, tau) / ccdf(d, tau)
            @test hazard(d, te, t) ≈ expected rtol=1e-10
        end
    end

    @testset "general UnivariateDistribution (Gamma)" begin
        # Test with Gamma distribution
        alpha = 2.0  # shape
        scale = 1.5
        d = Gamma(alpha, scale)
        te = 0.0

        for t in [0.1, 1.0, 2.0, 5.0]
            expected = pdf(d, t) / ccdf(d, t)
            @test hazard(d, te, t) ≈ expected rtol=1e-10
        end
    end

    @testset "Dual-parameterized hazard value part equals Float64 hazard" begin
        # Loosening the time signature to ::Real lets ForwardDiff.Dual
        # distribution parameters flow through hazard. The Dual's value part must
        # reproduce the plain Float64 hazard at the same parameters, since the
        # value part is exactly the undifferentiated computation.
        k, scale = 2.3, 3.1
        te, t = 1.0, 4.7

        function h(θ)
            d = Weibull(θ[1], θ[2])
            return hazard(d, te, t)
        end
        θ = [k, scale]
        dual_result = ForwardDiff.gradient(h, θ)  # forces Dual evaluation inside
        @test dual_result isa Vector{Float64}

        # Evaluate hazard with a Dual directly and compare its value part.
        d_dual = Weibull(ForwardDiff.Dual(k, 1.0, 0.0), ForwardDiff.Dual(scale, 0.0, 1.0))
        hd = hazard(d_dual, te, t)
        @test hd isa ForwardDiff.Dual
        @test ForwardDiff.value(hd) ≈ hazard(Weibull(k, scale), te, t) rtol=1e-12
    end
end
