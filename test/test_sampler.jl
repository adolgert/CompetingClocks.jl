using SafeTestsets


module MultiSamplerHelp
using CompetingClocks
using CompetingClocks: DirectCall, FirstToFire, MultiSampler, SamplerChoice, choose_sampler, jitter!
using Distributions: Exponential, UnivariateDistribution
using Random: AbstractRNG

struct ByDistribution <: SamplerChoice{Int64,Int64} end

function CompetingClocks.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::Exponential
)::Int64
    return 1
end
function CompetingClocks.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::UnivariateDistribution
)::Int64
    return 2
end

struct ByRate <: SamplerChoice{String,Int64} end

function CompetingClocks.choose_sampler(
    chooser::ByRate, clock::Int64, distribution::UnivariateDistribution
)::String
    if clock ≥ 50 && clock < 100
        return "fast"
    else
        return "slow"
    end
end

# A mock sampler that records the `when` argument it is queried with in `next`
# and returns a fixed (time, key). Used to verify time propagation in
# MultiSampler.next.
mutable struct WhenRecorder{K,T} <: CompetingClocks.SSA{K,T}
    recorded::T
    fire_time::T
    fire_key::K
end

function CompetingClocks.next(s::WhenRecorder{K,T}, when::T, rng::AbstractRNG) where {K,T}
    s.recorded = when
    return (s.fire_time, s.fire_key)
end
end


@safetestset multisampler_smoke = "MultiSampler smoke" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!, disable!, next, choose_sampler, reset!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential, Gamma

    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()
    rng = Xoshiro(90422342)
    enabled = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        @debug "Calling enable on $clock_id"
        if clock_id < 3
            enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0, rng)
        else
            enable!(sampler, clock_id, Gamma(propensity), 0.0, 0.0, rng)
        end
        push!(enabled, clock_id)
    end
    @test 1 ∈ keys(sampler.propagator[1])
    @test 2 ∈ keys(sampler.propagator[1])
    @test 3 ∈ keys(sampler.propagator[2])

    @test haskey(sampler, 1)
    @test !haskey(sampler, 1_000)
    @test !haskey(sampler, "1")

    when, which = next(sampler, 0.0, rng)
    @test which !== nothing
    disable!(sampler, which, when)
    delete!(enabled, which)
    @test enabled == Set(keys(sampler))
    todisable = pop!(enabled)
    disable!(sampler, todisable, when)
    @test enabled == Set(keys(sampler))
    enable!(sampler, 35, Exponential(), when, when, rng)
    push!(enabled, 35)
    when, which = next(sampler, when, rng)
    @test which ∈ enabled
    when, which = next(sampler, when, rng)
    @test which ∈ enabled
    reset!(sampler)
end


@safetestset multisampler_hierarchical = "MultiSampler hierarchical" begin
    # Let's make a hierarchical sampler that contains a hierarchical sampler.
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, DirectCall, MultiSampler, enable!
    using CompetingClocks: disable!, next, choose_sampler, reset!, enabled
    using ..MultiSamplerHelp: ByDistribution, ByRate
    using Distributions: Exponential, Gamma

    EventKey = Int64
    Time = Float64
    sampler = MultiSampler{Int64,EventKey,Time}(ByDistribution())
    sampler[1] = DirectCall{EventKey,Time}()
    sampler[2] = FirstToFire{EventKey,Time}()
    highest = MultiSampler{String,EventKey,Time}(ByRate())
    highest["fast"] = FirstToFire{EventKey,Time}()
    highest["slow"] = sampler

    rng = Xoshiro(90422342)
    enabled_set = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        @debug "Calling enable on $clock_id"
        if clock_id < 3
            enable!(highest, clock_id, Exponential(propensity), 0.0, 0.0, rng)
        else
            enable!(highest, clock_id, Gamma(propensity), 0.0, 0.0, rng)
        end
        push!(enabled_set, clock_id)
    end
    @test length(enabled(sampler)) == 5
    @test 3 in enabled(sampler)
    @test 149 ∉ enabled(sampler)
    enable!(highest, 53, Exponential(10.0), 0.0, 0.0, rng)
    push!(enabled_set, 53)
    enable!(highest, 57, Exponential(10.0), 0.0, 0.0, rng)
    push!(enabled_set, 57)
    @test 1 ∈ keys(highest.propagator["slow"].propagator[1])
    @test 2 ∈ keys(highest.propagator["slow"].propagator[1])
    @test 3 ∈ keys(highest.propagator["slow"].propagator[2])
    @test 53 ∈ keys(highest.propagator["fast"])
    when, which = next(highest, 0.0, rng)
    @test which !== nothing
    disable!(highest, which, when)
    delete!(enabled_set, which)
    @test enabled_set == Set(keys(highest))
    todisable = pop!(enabled_set)
    disable!(highest, todisable, when)
    @test enabled_set == Set(keys(highest))
    enable!(highest, 35, Exponential(), when, when, rng)
    push!(enabled_set, 35)
    when, which = next(highest, when, rng)
    @test which ∈ enabled_set
    when, which = next(highest, when, rng)
    @test which ∈ enabled_set
    reset!(highest)
end


@safetestset multisampler_clone = "MultiSampler clone" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!, clone
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential

    rng = Xoshiro(234567)
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)

    cloned = clone(sampler)
    @test length(cloned.propagator) == 2  # same structure
    @test haskey(cloned.propagator, 1)
    @test haskey(cloned.propagator, 2)
    @test length(cloned) == 0  # cloned subsamplers are empty
    @test length(sampler) == 2  # original unchanged
end


@safetestset multisampler_copy_clocks = "MultiSampler copy_clocks!" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!, copy_clocks!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential

    rng = Xoshiro(345678)

    src = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    src[1] = FirstToFire{Int64,Float64}()
    src[2] = FirstToFire{Int64,Float64}()
    enable!(src, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(src, 2, Exponential(2.0), 0.0, 0.0, rng)

    dst = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    dst[1] = FirstToFire{Int64,Float64}()
    dst[2] = FirstToFire{Int64,Float64}()

    @test length(src) == 2
    @test length(dst) == 0

    copy_clocks!(dst, src)
    # After copy, dst should have same state
    @test length(dst.chosen) == 2
end


@safetestset multisampler_jitter = "MultiSampler jitter!" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!, next
    using CompetingClocks: jitter!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential

    rng = Xoshiro(456789)
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)

    # Get initial next event
    t1, k1 = next(sampler, 0.0, rng)

    # Jitter should resample - times may change
    jitter!(sampler, 0.5, rng)

    t2, k2 = next(sampler, 0.5, rng)
    @test t2 >= 0.5
end


@safetestset multisampler_fire = "MultiSampler fire!" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!, next, fire!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential

    rng = Xoshiro(567890)
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)
    @test length(sampler) == 2

    when, which = next(sampler, 0.0, rng)
    fire!(sampler, which, when)
    @test length(sampler) == 1
end


@safetestset multisampler_getindex = "MultiSampler getindex" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential

    rng = Xoshiro(678901)
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)

    # getindex returns the stored firing time from FirstToFire
    firing_time = sampler[1]
    @test firing_time > 0.0
    @test firing_time < Inf
end


@safetestset multisampler_length = "MultiSampler length" begin
    using Random: Xoshiro
    using CompetingClocks: FirstToFire, MultiSampler, enable!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential, Gamma

    rng = Xoshiro(789012)
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()

    @test length(sampler) == 0

    enable!(sampler, 1, Exponential(1.0), 0.0, 0.0, rng)
    @test length(sampler) == 1

    enable!(sampler, 2, Exponential(2.0), 0.0, 0.0, rng)
    @test length(sampler) == 2

    # Gamma goes to sampler[2]
    enable!(sampler, 3, Gamma(1.0, 1.0), 0.0, 0.0, rng)
    @test length(sampler) == 3
end


@safetestset multisampler_choose_fallback = "MultiSampler choose_sampler fallback" begin
    using CompetingClocks: SamplerChoice, choose_sampler
    using Distributions: Exponential

    # A chooser that doesn't implement choose_sampler
    struct NoImplementation <: SamplerChoice{Symbol,Int64} end

    chooser = NoImplementation()
    @test_throws MissingException choose_sampler(chooser, 1, Exponential(1.0))
end


@safetestset multisampler_next_time_propagation = "MultiSampler next time propagation" begin
    using Random: Xoshiro
    using CompetingClocks: MultiSampler, next
    using ..MultiSamplerHelp: ByDistribution, WhenRecorder

    # Two sub-samplers, each recording the `when` it is queried with. They
    # return different fixed firing times, so the loop's `least_when` changes
    # between iterations. The bug reassigned `when`, causing the second
    # sub-sampler to be queried with the first's firing time.
    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    rec1 = WhenRecorder{Int64,Float64}(NaN, 3.0, 11)
    rec2 = WhenRecorder{Int64,Float64}(NaN, 7.0, 22)
    sampler[1] = rec1
    sampler[2] = rec2

    rng = Xoshiro(12345)
    current_time = 5.0
    when, which = next(sampler, current_time, rng)

    # Both sub-samplers must have been queried with the true current time.
    @test rec1.recorded == current_time
    @test rec2.recorded == current_time
    # The soonest firing time is returned.
    @test when == 3.0
    @test which == 11
end
