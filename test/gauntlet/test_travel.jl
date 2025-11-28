using Test
using CompetingClocks
using Random
using Distributions
using Graphs
using Base.Threads
using HypothesisTests

include("travel.jl")
using .TravelModel
include("generate_data.jl")
include("mark_calibration.jl")
include("running_score.jl")
include("experiments.jl")



# Helper to create a vector of RNGs for threaded operations
function make_rng_vector(seed::Integer)
    rng = Vector{Xoshiro}(undef, max(Threads.maxthreadid(), 1))
    rng[1] = Xoshiro(seed)
    for i in 2:length(rng)
        rng[i] = Xoshiro(seed + i)
    end
    return rng
end


@testset "Travel Model Smoke Tests" begin
    @testset "Graph Construction" begin
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

    @testset "TravelModel Construction" begin
        rng = Xoshiro(12345)

        config = TravelConfig(
            TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model = Travel(3, config, rng)
        @test model isa Travel
        @test nv(model.g) == 3
        @test !isempty(model.rates)
        @test model.remember == TravelMemory.forget

        config_remember = TravelConfig(
            TravelMemory.remember, TravelGraph.path, TravelRateDist.exponential,
            TravelRateCount.destination, TravelRateDelay.none
            )
        model_remember = Travel(3, config_remember, rng)
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

        rec = VectorRecord()
        travel_run(5, sampler, model, (x,y) -> nothing, rec, rng)
        commands = rec.commands
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
        rec = VectorRecord()
        travel_run(5, sampler, model, (x,y) -> nothing, rec, rng)

        sampler2 = FirstReaction{Int,Float64}()
        rng2 = Xoshiro(111)
        replay_commands(rec.commands, sampler2, rng2)
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
        rec = VectorRecord()
        travel_run(5, sampler, model, (x,y) -> nothing, rec, rng)
        commands = rec.commands

        rng_vec = make_rng_vector(12345)
        samplers, when = parallel_replay(commands, sampler, 4, rng_vec)
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
        rec = VectorRecord()
        travel_run(5, sampler, model, (x,y) -> nothing, rec, rng)
        commands = rec.commands
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
        rec = VectorRecord()
        travel_run(5, sampler, model, (x,y) -> nothing, rec, rng)
        commands = rec.commands

        rng_vec = make_rng_vector(12345)
        samplers, when = parallel_replay(commands, sampler, 4, rng_vec)
        rng_data = make_rng_vector(999)
        data = sample_samplers(samplers, when, rng_data)

        @test data isa Vector
        @test length(data) == 4
        @test all(d -> d isa Tuple{Int,Float64}, data)
        @test all(d -> d[2] > when, data)
    end
end


@testset "Mark Calibration Smoke Tests" begin
    @testset "Mark Calibration Brier" begin
        # Create simple test case using DistributionState
        dist1 = Exponential(1.0)
        dist2 = Exponential(2.0)
        distributions = Dict(
            1 => DistributionState(dist1, 0.0),
            2 => DistributionState(dist2, 0.0)
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
        rec = VectorRecord()
        travel_run(10, sampler, model, (x,y) -> nothing, rec, rng)
        commands = rec.commands

        rng_vec = make_rng_vector(12345)
        samplers, when = parallel_replay(commands, sampler, 4, rng_vec)
        rng_data = make_rng_vector(999)
        draws = sample_samplers(samplers, when, rng_data)

        final_dist = final_enabled_distributions(commands)

        # Test that it runs and returns a number
        total = mark_calibration_conditional_time(draws, final_dist)
        @test total isa Float64
        @test total >= 0.0
        @test isfinite(total)
    end
end


@testset "Experiments Integration Tests" begin
    @testset "run_experiments faster mode" begin
        @test run_experiments(true)
    end
end
