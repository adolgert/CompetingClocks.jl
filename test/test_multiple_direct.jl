using SafeTestsets

module MultipleKeyedHelp
using CompetingClocks
using Distributions: Exponential, UnivariateDistribution

struct ByName <: SamplerChoice{String,Symbol} end
export ByName

function CompetingClocks.choose_sampler(
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
    using CompetingClocks
    using ..MultipleKeyedHelp

    rng = Xoshiro(432234)
    SamplerKey = Symbol # The samplers
    K = String  # The clocks
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    cumulant_scan = CompetingClocks.CumSumPrefixSearch{Time}()
    keep_prefix_tree = CompetingClocks.KeyedKeepPrefixSearch{K,typeof(cumulant_scan)}(cumulant_scan)
    md[:fast] = keep_prefix_tree

    when = 0.0
    enable!(md, "moveleft", Exponential(0.35), 0.0, when, rng)
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

    @test length(enabled(md)) == 2
    # Disabling from the removal tree will remove it.
    disable!(md, "moveleft", when)
    @test length(keyed_prefix_tree) == 0
    @test enabled(md) == Set(["fastlighton"])

    # Disabling from the keep tree won't remove it.
    disable!(md, "fastlighton", when)
    @test length(keep_prefix_tree) == 1
    reset!(md)
end


@safetestset multiple_direct_clone = "MultipleDirect clone" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using ..MultipleKeyedHelp

    rng = Xoshiro(234567)
    SamplerKey = Symbol
    K = String
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    enable!(md, "moveleft", Exponential(0.35), 0.0, 0.0, rng)

    cloned = clone(md)
    # Cloned has empty scan/totals/chosen/scanmap vectors
    @test length(cloned.scan) == 0
    @test length(cloned.chosen) == 0
    @test length(enabled(md)) == 1  # original unchanged
end


@safetestset multiple_direct_jitter = "MultipleDirect jitter!" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using CompetingClocks: jitter!
    using ..MultipleKeyedHelp

    rng = Xoshiro(345678)
    SamplerKey = Symbol
    K = String
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    enable!(md, "moveleft", Exponential(0.35), 0.0, 0.0, rng)

    # jitter! does nothing for MultipleDirect (returns nothing)
    result = jitter!(md, 0.5, rng)
    @test result === nothing
end


@safetestset multiple_direct_copy = "MultipleDirect copy_clocks!" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using ..MultipleKeyedHelp

    rng = Xoshiro(456789)
    SamplerKey = Symbol
    K = String
    Time = Float64

    # Create source
    src = MultipleDirect{SamplerKey,K,Time}(ByName())
    src_prefix = CompetingClocks.BinaryTreePrefixSearch{Time}()
    src_keyed = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(src_prefix)}(src_prefix)
    src[:slow] = src_keyed
    src_fast_prefix = CompetingClocks.CumSumPrefixSearch{Time}()
    src_fast_keyed = CompetingClocks.KeyedKeepPrefixSearch{K,typeof(src_fast_prefix)}(src_fast_prefix)
    src[:fast] = src_fast_keyed

    enable!(src, "moveleft", Exponential(0.35), 0.0, 0.0, rng)
    enable!(src, "fastaction", Exponential(1.0), 0.0, 0.0, rng)

    # Create destination with same structure
    dst = MultipleDirect{SamplerKey,K,Time}(ByName())
    dst_prefix = CompetingClocks.BinaryTreePrefixSearch{Time}()
    dst_keyed = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(dst_prefix)}(dst_prefix)
    dst[:slow] = dst_keyed
    dst_fast_prefix = CompetingClocks.CumSumPrefixSearch{Time}()
    dst_fast_keyed = CompetingClocks.KeyedKeepPrefixSearch{K,typeof(dst_fast_prefix)}(dst_fast_prefix)
    dst[:fast] = dst_fast_keyed

    @test length(enabled(src)) == 2
    @test length(enabled(dst)) == 0

    copy_clocks!(dst, src)
    @test length(enabled(dst)) == 2
end


@safetestset multiple_direct_reenable = "MultipleDirect re-enable clock" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using ..MultipleKeyedHelp

    rng = Xoshiro(567890)
    SamplerKey = Symbol
    K = String
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    # Enable a clock
    enable!(md, "moveleft", Exponential(0.35), 0.0, 0.0, rng)
    @test length(enabled(md)) == 1

    # Re-enable the same clock with different rate (exercises line 95)
    enable!(md, "moveleft", Exponential(0.5), 0.0, 0.1, rng)
    @test length(enabled(md)) == 1  # Still just one clock
end


@safetestset multiple_direct_fire = "MultipleDirect fire!" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using ..MultipleKeyedHelp

    rng = Xoshiro(678901)
    SamplerKey = Symbol
    K = String
    Time = Float64

    md = MultipleDirect{SamplerKey,K,Time}(ByName())
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree
    fast_prefix = CompetingClocks.CumSumPrefixSearch{Time}()
    fast_keyed = CompetingClocks.KeyedKeepPrefixSearch{K,typeof(fast_prefix)}(fast_prefix)
    md[:fast] = fast_keyed

    enable!(md, "moveleft", Exponential(0.35), 0.0, 0.0, rng)
    enable!(md, "fastaction", Exponential(1.0), 0.0, 0.0, rng)
    @test length(enabled(md)) == 2

    # Get next event and fire it
    when, which = next(md, 0.0, rng)
    fire!(md, which, when)

    # One clock should be removed (the slow one if it was chosen)
    # or zeroed out (the fast one if it was chosen from keep tree)
    @test length(enabled(md)) == 1
end


@safetestset multiple_direct_trajectory = "MultipleDirect trajectory likelihood" begin
    using Distributions: Exponential
    using Random: Xoshiro
    using CompetingClocks
    using CompetingClocks: steploglikelihood, pathloglikelihood
    using ..MultipleKeyedHelp

    rng = Xoshiro(890123)
    SamplerKey = Symbol
    K = String
    Time = Float64

    # Create with trajectory=true
    md = MultipleDirect{SamplerKey,K,Time}(ByName(), trajectory=true)
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Time}()
    keyed_prefix_tree = CompetingClocks.KeyedRemovalPrefixSearch{K,typeof(prefix_tree)}(prefix_tree)
    md[:slow] = keyed_prefix_tree

    enable!(md, "moveleft", Exponential(0.5), 0.0, 0.0, rng)  # rate = 2.0
    enable!(md, "moveright", Exponential(1.0), 0.0, 0.0, rng)  # rate = 1.0

    # Test steploglikelihood
    now = 0.0
    when = 0.1
    which = "moveleft"
    ll = steploglikelihood(md, now, when, which)
    # log(2.0) - 3.0 * 0.1 = 0.693... - 0.3 ≈ 0.393
    @test ll ≈ log(2.0) - 3.0 * 0.1

    # Fire the event to update md.now and md.log_likelihood
    fire!(md, which, when)

    # Test pathloglikelihood
    path_ll = pathloglikelihood(md, 0.2)
    @test isfinite(path_ll)
end
