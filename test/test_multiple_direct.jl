using SafeTestsets

module MultipleKeyedHelp
    using Fleck
    using Distributions: Exponential, UnivariateDistribution

    struct ByName <: SamplerChoice{String,Symbol} end
    export ByName

    function Fleck.choose_sampler(
        chooser::ByName, clock::String, distribution::UnivariateDistribution
        )::Symbol
        if startswith(clock, "fast")
            return :fast
        else
            return :slow
        end
    end
end


@safetestset multiple_direct_smoke = "MultipleDirect smoke" begin
    using Distributions
    using Random
    using Fleck
    using ..MultipleKeyedHelp

    rng = Xoshiro(432234)
    SamplerKey = Symbol # The samplers
    K = String  # The clocks
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = Fleck.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = Fleck.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    cumulant_scan = Fleck.CumSumPrefixSearch{Time}()
    keep_prefix_tree = Fleck.KeyedKeepPrefixSearch{K,typeof(cumulant_scan)}(cumulant_scan)
    md[:fast] = keep_prefix_tree

    when = 0.0
    enable!(md, "moveleft", Exponential(.35), 0.0, when, rng)
    @test length(keyed_prefix_tree) == 1
    @test length(keep_prefix_tree) == 0

    enable!(md, "fastlighton", Exponential(1.0), 0.0, when, rng)
    @test length(keyed_prefix_tree) == 1
    @test length(keep_prefix_tree) == 1

    seen = Set{String}()
    for sample_idx in 1:100
        when, which = next(md, 0.0, rng)
        push!(seen, which)
    end
    @test "moveleft" ∈ seen
    @test "fastlighton" ∈ seen

    # Disabling from the removal tree will remove it.
    disable!(md, "moveleft", when)
    @test length(keyed_prefix_tree) == 0

    # Disabling from the keep tree won't remove it.
    disable!(md, "fastlighton", when)
    @test length(keep_prefix_tree) == 1
end
