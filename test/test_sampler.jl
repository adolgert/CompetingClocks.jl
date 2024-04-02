using SafeTestsets

@safetestset singlesampler_smoke = "SingleSampler smoke" begin
    using Random: Xoshiro
    using Fleck: FirstToFire, SingleSampler, enable!, disable!, sample!
    using Distributions: Exponential

    sampler = SingleSampler(FirstToFire{Int64,Float64}())
    rng = Xoshiro(90422342)
    enabled = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity), 0.0, rng)
        push!(enabled, clock_id)
    end
    when, which = sample!(sampler, rng)
    delete!(enabled, when)
    todisable = collect(enabled)[1]
    disable!(sampler, todisable)
    enable!(sampler, 35, Exponential(), when, rng)
    when, which = sample!(sampler, rng)
end


module MultiSamplerHelp
    using Fleck
    using Distributions: Exponential, UnivariateDistribution

    struct ByDistribution <: SamplerChoice{Int64,Int64} end

    function Fleck.choose_sampler(
        chooser::ByDistribution, clock::Int64, distribution::Exponential
        )::Int64
        return 1
    end
    function Fleck.choose_sampler(
        chooser::ByDistribution, clock::Int64, distribution::UnivariateDistribution
        )::Int64
        return 2
    end
end


@safetestset multisampler_smoke = "MultiSampler smoke" begin
    using Random: Xoshiro
    using Fleck: FirstToFire, MultiSampler, enable!, disable!, sample!, choose_sampler, reset!
    using ..MultiSamplerHelp: ByDistribution
    using Distributions: Exponential, Gamma

    sampler = MultiSampler{Int64,Int64,Float64}(ByDistribution())
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()
    rng = Xoshiro(90422342)
    enabled = Set{Int64}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        if clock_id < 3
            enable!(sampler, clock_id, Exponential(propensity), 0.0, rng)
        else
            enable!(sampler, clock_id, Gamma(propensity), 0.0, rng)
        end
        push!(enabled, clock_id)
    end
    when, which = sample!(sampler, rng)
    delete!(enabled, when)
    todisable = collect(enabled)[1]
    disable!(sampler, todisable)
    enable!(sampler, 35, Exponential(), when, rng)
    when, which = sample!(sampler, rng)
    reset!(sampler)
end
