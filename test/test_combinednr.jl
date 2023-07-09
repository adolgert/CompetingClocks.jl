using SafeTestsets


@safetestset CombinedNextReactionSmoke = "combinednext reaction does basic things" begin
    using Distributions
    using Random
    using Fleck: CombinedNextReaction, next, enable!, disable!

    rng = MersenneTwister(349827)
    for i in 1:100
        sampler = CombinedNextReaction{String}()
        @test next(sampler, 3.0, rng)[2] === nothing
        enable!(sampler, "walk home", Exponential(1.5), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] == "walk home"
        enable!(sampler, "run", Gamma(1, 3), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        enable!(sampler, "walk to sandwich shop", Weibull(2, 1), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run", "walk to sandwich shop"]
        disable!(sampler, "walk to sandwich shop", 1.7)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
    end
end
