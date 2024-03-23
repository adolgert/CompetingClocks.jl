using SafeTestsets


@safetestset firsttofire_smoke = "FirstToFire smoke" begin
    using Fleck: FirstToFire, enable!, disable!, next
    using Random: Xoshiro
    using Distributions: Exponential

    sampler = FirstToFire{Int64,Float64}()
    rng = Xoshiro(90422342)
    enabled = Set{Int64}()
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


@safetestset firsttofire_insertion = "FirstToFire known insertions" begin
    # Let's put entries into its heap in order to ensure the heap does the right thing.
    # This is a clear-box test. It depends heavily on implementation.
    using Fleck: FirstToFire, enable!, disable!, next
    using Fleck
    using Random: Xoshiro

    propagator = FirstToFire{Int64,Float64}()
    
    @test length(propagator) == 0
    @test length(keys(propagator)) == 0
    @test_throws KeyError propagator[1]
    @test keytype(propagator) <: Int64
    
    for (clock, when_fire) in [(1, 7.9), (2, 12.3), (3, 3.7), (4, 0.00013), (5, 0.2)]
        heap_handle = push!(
            propagator.firing_queue,
            Fleck.OrderedSample{Int64,Float64}(clock, when_fire)
            )
        propagator.transition_entry[clock] = heap_handle
    end

    @test length(propagator) == 5
    @test length(keys(propagator)) == 5
    @test propagator[1] == 7.9

    rng = Xoshiro(39472)
    for (key, fire_time) in [(4, 0.00013),(5, 0.2), (3, 3.7), (1, 7.9), (2, 12.3)]
        (t1, k1) = next(propagator, 0.0, rng)
        @test k1 == key
        @test abs(t1 - fire_time) < 0.00001
        disable!(propagator, k1, 0.0)
    end
    (t2, k2) = next(propagator, 0.0, rng)
    @test k2 === nothing
end
