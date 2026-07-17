using CompetingClocks

@safetestset lefttrunc_smoke = "Left-truncated exponential rand" begin
    using Distributions
    using Random

    @testset "lower-only truncation leaves upper === nothing" begin
        rng = Xoshiro(979797)
        d = truncated(Exponential(1.5); lower=2.0)
        for _ in 1:100
            draw = rand(rng, d)
            @test draw >= 2.0
        end
    end

    @testset "explicit infinite upper bound" begin
        rng = Xoshiro(24601)
        d = truncated(Exponential(1.5), 2.0, Inf)
        for _ in 1:100
            draw = rand(rng, d)
            @test draw >= 2.0
        end
    end

    @testset "memoryless shift matches untruncated draw" begin
        d = truncated(Exponential(3.0); lower=5.0)
        @test rand(Xoshiro(42), d) == 5.0 + rand(Xoshiro(42), Exponential(3.0))
    end
end
