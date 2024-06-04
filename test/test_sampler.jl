using SafeTestsets


module MultiSamplerHelp
    using CompetingClocks
    using Distributions: Exponential, UnivariateDistribution

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
    using CompetingClocks: disable!, next, choose_sampler, reset!
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
    enabled = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        @debug "Calling enable on $clock_id"
        if clock_id < 3
            enable!(highest, clock_id, Exponential(propensity), 0.0, 0.0, rng)
        else
            enable!(highest, clock_id, Gamma(propensity), 0.0, 0.0, rng)
        end
        push!(enabled, clock_id)
    end
    enable!(highest, 53, Exponential(10.0), 0.0, 0.0, rng)
    push!(enabled, 53)
    enable!(highest, 57, Exponential(10.0), 0.0, 0.0, rng)
    push!(enabled, 57)
    @test 1 ∈ keys(highest.propagator["slow"].propagator[1])
    @test 2 ∈ keys(highest.propagator["slow"].propagator[1])
    @test 3 ∈ keys(highest.propagator["slow"].propagator[2])
    @test 53 ∈ keys(highest.propagator["fast"])
    when, which = next(highest, 0.0, rng)
    @test which !== nothing
    disable!(highest, which, when)
    delete!(enabled, which)
    @test enabled == Set(keys(highest))
    todisable = pop!(enabled)
    disable!(highest, todisable, when)
    @test enabled == Set(keys(highest))
    enable!(highest, 35, Exponential(), when, when, rng)
    push!(enabled, 35)
    when, which = next(highest, when, rng)
    @test which ∈ enabled
    when, which = next(highest, when, rng)
    @test which ∈ enabled
    reset!(highest)
end
