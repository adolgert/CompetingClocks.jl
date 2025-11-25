using Test
using CompetingClocks
using Random
using Distributions
using Graphs
using Base.Threads

include("travel.jl")
include("generate_data.jl")
include("single_clock.jl")

using .TravelModel


@testset "Travel Model Smoke Tests" begin
    @testset "Graph Construction" begin
        rng = Xoshiro(12345)

        # Test each graph type
        g1 = travel_make_graph(TravelGraph.path, 5)
        @test g1 isa SimpleGraph
        @test nv(g1) == 5

        g2 = travel_make_graph(TravelGraph.cycle, 5)
        @test g2 isa SimpleGraph
        @test nv(g2) == 5

        g3 = travel_make_graph(TravelGraph.complete, 5)
        @test g3 isa SimpleGraph
        @test nv(g3) == 5

        g4 = travel_make_graph(TravelGraph.clique, 5)
        @test g4 isa SimpleGraph
    end

    @testset "Rate Generation" begin
        rng = Xoshiro(12345)

        rate_tuple = travel_make_rate(1, rng)
        @test rate_tuple isa Tuple
        @test length(rate_tuple) == 2
        @test rate_tuple[1] isa UnivariateDistribution
        @test rate_tuple[2] isa Float64

        rates = travel_rates_exponential(3, rng)
        @test rates isa Dict
        @test !isempty(rates)
    end

    @testset "TravelModel Construction" begin
        rng = Xoshiro(12345)

        config = TravelConfig(
            TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, TravelGraph.path, TravelMemory.forget, rng)
        @test model isa Travel
        @test nv(model.g) == 3
        @test !isempty(model.rates)
        @test model.remember == TravelMemory.forget

        model_remember = Travel(3, TravelGraph.path, TravelMemory.remember, rng)
        @test model_remember.remember == TravelMemory.remember
    end

    @testset "State Initialization and Enabled" begin
        rng = Xoshiro(12345)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.cycle, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, config, rng)

        state = travel_init_state(model, rng)
        @test state isa Int
        @test 1 <= state <= 3

        enabled = travel_enabled(model, state)
        @test enabled isa Set
    end

    @testset "Simulation Run" begin
        rng = Xoshiro(12345)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.cycle, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, config, rng)
        sampler = FirstReaction{Int,Float64}()

        commands = travel_run(5, sampler, model, rng)
        @test commands isa Vector
        @test !isempty(commands)
        @test all(cmd -> cmd isa Tuple, commands)
        @test all(cmd -> cmd[1] in (:enable, :disable, :fire), commands)
    end

    @testset "Complete Workflow" begin
        # This tests travel_make_run()
        result = travel_make_run()
        # Should complete without error
        @test true
    end
end


@testset "Generate Data Smoke Tests" begin
    @testset "Replay Commands" begin
        rng = Xoshiro(98327423)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(2, config, rng)
        sampler = FirstReaction{Int,Float64}()
        commands = travel_run(5, sampler, model, rng)

        sampler2 = FirstReaction{Int,Float64}()
        rng2 = Xoshiro(111)
        replay_commands(commands, sampler2, rng2)
        # Should complete without error
        @test true
    end

    @testset "Parallel Replay" begin
        rng = Xoshiro(98327423)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(2, config, rng)
        sampler = FirstReaction{Int,Float64}()
        commands = travel_run(5, sampler, model, rng)

        samplers, when = parallel_replay(commands, 4, 12345)
        @test samplers isa Vector
        @test length(samplers) == 4
        @test all(s -> s isa FirstReaction, samplers)
        @test when isa Float64
        @test when > 0
    end

    @testset "Final Enabled Distributions" begin
        rng = Xoshiro(98327423)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.cycle, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, config, rng)
        sampler = FirstReaction{Int,Float64}()
        commands = travel_run(5, sampler, model, rng)
        @test !isempty(commands)

        final_dist = final_enabled_distributions(commands)
        @test final_dist isa Dict
        # Should have at least one enabled distribution
        @test !isempty(final_dist)
        for (k, v) in final_dist
            @test k isa Int
            @test v isa DistributionState
            @test v.d isa UnivariateDistribution
            @test v.enabling_time isa Float64
        end
    end

    @testset "Sample Samplers" begin
        rng = Xoshiro(98327423)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(2, config, rng)
        sampler = FirstReaction{Int,Float64}()
        commands = travel_run(5, sampler, model, rng)

        samplers, when = parallel_replay(commands, 4, 12345)
        rng_data = Xoshiro(999)
        data = sample_samplers(samplers, when, rng_data)

        @test data isa Vector
        @test length(data) == 4
        @test all(d -> d isa Tuple{Int,Float64}, data)
        @test all(d -> d[2] > when, data)
    end
end


@testset "Single Clock Smoke Tests" begin
    @testset "Mark Calibration Brier" begin
        # Create simple test case
        dist1 = Exponential(1.0)
        dist2 = Exponential(2.0)
        distributions = Dict(
            1 => (dist1, 0.0),
            2 => (dist2, 0.0)
        )

        # Test that it runs and returns a number
        brier = mark_calibration_brier(distributions, 1, 1.0)
        @test brier isa Float64
        @test brier >= 0.0
        @test isfinite(brier)
    end

    @testset "Mark Calibration Conditional Time" begin
        # Generate some test data
        rng = Xoshiro(98327423)
        config = TravelConfig(
            TravelMemory.forget, TravelGraph.cycle, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, config, rng)
        sampler = FirstReaction{Int,Float64}()
        commands = travel_run(10, sampler, model, rng)

        samplers, when = parallel_replay(commands, 4, 12345)
        rng_data = Xoshiro(999)
        draws = sample_samplers(samplers, when, rng_data)

        final_dist = final_enabled_distributions(commands)
        # Convert to expected format
        distributions = Dict(k => (v.d, v.enabling_time) for (k, v) in final_dist)

        # Test that it runs and returns a number
        total = mark_calibration_conditional_time(draws, distributions)
        @test total isa Float64
        @test total >= 0.0
        @test isfinite(total)
    end
end
