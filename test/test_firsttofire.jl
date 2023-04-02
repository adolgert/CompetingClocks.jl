using SafeTestsets


@safetestset firsttofire_smoke = "FirstToFire smoke" begin
    using Fleck: FirstToFire, enable!, disable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    sampler = FirstToFire{Int}()
    rng = Xoshiro(90422342)
    enabled = Set{Int}()
    for (clock_id, propensity) in enumerate([0.3, 0.2, 0.7, 0.001, 0.25])
        enable!(sampler, clock_id, Exponential(propensity), 0.0, 0.0, rng)
        push!(enabled, clock_id)
    end
    when, which = next(sampler, 0.0, rng)
    disable!(sampler, which, when)
    delete!(enabled, which)
    when, which = next(sampler, when, rng)
    @test when > 0.0
    @test 1 <= which
    @test which <= 5
    @test which ∈ enabled
end
