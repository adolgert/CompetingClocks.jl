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


@safetestset multisampler_smoke = "MultiSampler smoke" begin
    using Random: Xoshiro
    using Fleck: FirstToFire, MultiSampler, enable!, disable!, sample!
    using Distributions: Exponential

    chooser = clock_id -> begin
        if clock_id < 4
            return 1
        else
            return 2
        end
    end

    sampler = MultiSampler{Int64,Int64,Float64}(chooser)
    sampler[1] = FirstToFire{Int64,Float64}()
    sampler[2] = FirstToFire{Int64,Float64}()
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
