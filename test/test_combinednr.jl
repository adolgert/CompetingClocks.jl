using SafeTestsets


@safetestset CombinedNextReactionSmoke = "CombinedNextReaction reaction does basic things" begin
    using Distributions
    using Random
    using CompetingClocks: CombinedNextReaction, next, enable!, disable!, reset!

    rng = MersenneTwister(349827)
    for i in 1:100
        sampler = CombinedNextReaction{String,Float64}()
        @test next(sampler, 3.0, rng)[2] === nothing
        enable!(sampler, "walk home", Exponential(1.5), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] == "walk home"
        enable!(sampler, "run", Gamma(1, 3), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        enable!(sampler, "walk to sandwich shop", Weibull(2, 1), 0.0, 0.0, rng)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run", "walk to sandwich shop"]
        disable!(sampler, "walk to sandwich shop", 1.7)
        @test next(sampler, 3.0, rng)[2] ∈ ["walk home", "run"]
        reset!(sampler)
    end
end

@safetestset CombinedNextReaction_interface = "CombinedNextReaction basic interface" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    sampler = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[1]
    @test keytype(sampler) <: Int64

    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[1] == 7.9

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "1")

    disable!(sampler, 1, 0.0)

    @test_throws BoundsError sampler[1]
    @test sampler[2] == 12.3

end


@safetestset CombinedNextReaction_pair_keys = "CombinedNextReaction pair keys" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    # Make a key that looks like those in Gen.jl.
    KeyType = Pair{Symbol,Int64}
    sampler = CombinedNextReaction{KeyType,Float64}()
    rng = Xoshiro(123)

    @test length(sampler) == 0
    @test length(keys(sampler)) == 0
    @test_throws KeyError sampler[:S=>1]
    @test keytype(sampler) <: KeyType

    for (clock, when_fire) in [
        (:S => 1, 7.9), (:S => 2, 12.3), (:I => 3, 3.7), (:I => 4, 0.00013), (:S => 5, 0.2)
    ]
        enable!(sampler, clock, Dirac(when_fire), 0.0, 0.0, rng)
    end

    @test length(sampler) == 5
    @test length(keys(sampler)) == 5
    @test sampler[:S=>1] == 7.9

    @test haskey(sampler, :S => 1)
    @test !haskey(sampler, :I => 17)

    disable!(sampler, :S => 1, 0.0)

    @test_throws BoundsError sampler[:S=>1]
    @test sampler[:S=>2] == 12.3

end

@safetestset CombinedNextReaction_copy = "CombinedNextReaction copy" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    dst = clone(src)
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0, rng)
    enable!(src, 38, Exponential(), 0.0, 0.0, rng)
    enable!(dst, 29, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 2
    @test length(dst) == 1
    copy_clocks!(dst, src)
    @test length(src) == 2
    @test length(dst) == 2
    # Changing src doesn't change dst.
    enable!(src, 48, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 2
    # Changing dst doesn't change src.
    enable!(dst, 49, Exponential(), 0.0, 0.0, rng)
    @test length(src) == 3
    @test length(dst) == 3
end


@safetestset CombinedNextReaction_set = "CombinedNextReaction set" begin
    using CompetingClocks
    using Distributions
    using Random: Xoshiro

    src = CombinedNextReaction{Int64,Float64}()
    rng = Xoshiro(123)

    enable!(src, 37, Exponential(), 0.0, 0.0, rng)
    enable!(src, 38, Exponential(), 0.0, 0.0, rng)
    enable!(src, 29, Exponential(), 0.0, 0.0, rng)
    disable!(src, 37, 0.1)

    enabled_set = enabled(src)
    @test 38 in enabled_set
    @test 37 ∉ enabled_set
    @test length(enabled_set) == 2
    @test eltype(enabled_set) == Int64
    b = Set([x for x in enabled_set])
    @assert b == Set([38, 29])
end