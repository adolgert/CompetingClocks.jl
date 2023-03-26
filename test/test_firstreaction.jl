using SafeTestsets

module DirectFixture

end


@safetestset first_reaction_smoke = "FirstReaction smoke" begin
    using Fleck: FirstReaction, enable!, disable!, next
    using Random: MersenneTwister
    using Distributions: Exponential, Gamms

    rng = MersenneTwister(90422342)
    min_when = -1.0
    for i in 1:100
        fr = FirstReaction{Int}()
        enable!(fr, 1, Exponential(1.7), 0.0, 0.0, rng)
        enable!(fr, 2, Gamma(9, 0.5), 0.0, 0.0, rng)
        enable!(fr, 3, Gamma(2, 2.0), 0.0, 0.0, rng)
        curtime = 0.5
        disable!(fr, 2, curtime)
        when, which = next(md, curtime, rng)
        @test which ∈ [1, 3]
        @test when > curtime
        min_when = min(min_when, when)
    end
    @test min_when < curtime + 0.1
end
