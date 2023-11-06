using Fleck
using SafeTestsets

module CRNHelper
    export FakeSampler, next, enable!, disable!
    using Random
    using Fleck
    using Distributions
    # The first set of tests of commonrandom will construct a fake sampler.
    # That fake sampler will know whether it gets a unique rng starting point.
    # It will sometimes make a draw and sometimes not.

    # Assumptions: Clock IDs are integers. Don't care about times.
    struct FakeSampler
        draws::Vector{Tuple{Float64, Int64}}
        FakeSampler() = new(Vector{Tuple{Float64, Int64}}())
    end


    function Fleck.next(cr::FakeSampler, when::Float64, rng::AbstractRNG) where {Sampler}
        ((time, clock_id), element) = findmin(cr.draws)
        deleteat!(cr.draws, element)
        return (time, clock_id)
    end


    function Fleck.enable!(
        cr::FakeSampler, clock::T, distribution::UnivariateDistribution,
        te::Float64, when::Float64, rng::AbstractRNG) where {Sampler, T}
        push!(cr.draws, (when + rand(rng, distribution), clock))
    end


    function Fleck.disable!(cr::FakeSampler, clock::T, when::Float64) where {Sampler, T}
        toremove = findnext(x->x[2]==clock, cr.draws, 1)
        if toremove !== nothing
            deleteat!(cr.draws, toremove)
        end
    end
end


@safetestset crn_fake_sampler = "Common random numbers FakeSampler" begin
    using ..CRNHelper
    using Random: Xoshiro
    using Distributions

    rng = Xoshiro(324234)
    sampler = FakeSampler()
    enabled = Set(1:5)
    for startup in enabled
       enable!(sampler, startup, Exponential(), 0.0, 0.0, rng)
    end
    @test length(sampler.draws) == 5
    for out in 1:5
        (when, which) = next(sampler, 0.0, rng)
        pop!(enabled, which)
    end
    @test length(sampler.draws) == 0
end


@safetestset crn_recording = "Common random numbers recording" begin
    using Fleck
    using Random
    using ..CRNHelper
    using Distributions

    rng = Xoshiro(2934723)
    sampler = FakeSampler()
    record_sampler = CommonRandomRecorder(sampler, Int, Xoshiro)

    enabled = Set(1:5)
    for startup in enabled
       enable!(record_sampler, startup, Exponential(), 0.0, 0.0, rng)
    end
    @test length(sampler.draws) == 5
    for out in 1:5
        (when, which) = next(record_sampler, 0.0, rng)
        pop!(enabled, which)
    end
    @test length(sampler.draws) == 0
end


@safetestset crn_replay = "Common random numbers replay" begin
    using Fleck
    using Random
    using ..CRNHelper
    using Distributions

    rng = Xoshiro(2934723)
    sampler = FakeSampler()
    record_sampler = CommonRandomRecorder(sampler, Int, Xoshiro)

    enabled = Set(1:5)
    for startup in enabled
       enable!(record_sampler, startup, Exponential(), 0.0, 0.0, rng)
    end
    @test length(sampler.draws) == 5
    for out in 1:5
        (when, which) = next(record_sampler, 0.0, rng)
        pop!(enabled, which)
    end
    @test length(sampler.draws) == 0
end
