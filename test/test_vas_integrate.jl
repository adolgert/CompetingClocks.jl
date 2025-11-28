using SafeTestsets
using Test


@safetestset vas_sample_direct = "VAS samples with direct" begin
using Random: MersenneTwister
using CompetingClocks: DirectCall
using ..VectorAddition
using ..SampleVAS: sample_sir

    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = DirectCall{Int,Float64}()

    starting = zeros(Int, 3 * cnt)
    starting[2:cnt] .= 1
    starting[cnt + 1] = 1  # Start with one infected.
    initial_state = vas_initial(vas, starting)

    simulation = VectorAdditionFSM(vas, initial_state, sampler, rng)
    next_transition = -1
    step_cnt = 0
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
        step_cnt += 1
    end

    # Extract final state: S is indices 1:cnt, I is cnt+1:2cnt, R is 2cnt+1:3cnt
    state = simulation.state.state
    final_S = sum(state[1:cnt])
    final_I = sum(state[cnt+1:2*cnt])
    final_R = sum(state[2*cnt+1:3*cnt])

    # Conservation: total population unchanged
    @test final_S + final_I + final_R == cnt
    # Simulation ended because no more events (all infected recovered)
    @test final_I == 0
    # At least some infections occurred
    @test final_R >= 1
    # Step count = infections + recoveries = (cnt - final_S) + final_R
    @test step_cnt == (cnt - final_S) + final_R
end


@safetestset vas_sample_fr = "VAS samples with first reaction" begin
    using Random: MersenneTwister
    using CompetingClocks: FirstReaction
    using ..VectorAddition
    using ..SampleVAS: sample_sir

    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = FirstReaction{Int,Float64}()

    starting = zeros(Int, 3 * cnt)
    starting[2:cnt] .= 1
    starting[cnt + 1] = 1  # Start with one infected.
    initial_state = vas_initial(vas, starting)

    simulation = VectorAdditionFSM(vas, initial_state, sampler, rng)
    next_transition = -1
    step_cnt = 0
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
        step_cnt += 1
    end

    # Extract final state
    state = simulation.state.state
    final_S = sum(state[1:cnt])
    final_I = sum(state[cnt+1:2*cnt])
    final_R = sum(state[2*cnt+1:3*cnt])

    @test final_S + final_I + final_R == cnt
    @test final_I == 0
    @test final_R >= 1
    @test step_cnt == (cnt - final_S) + final_R
end



@safetestset vas_sample_ftf = "VAS samples with first to fire" begin
    using Random: MersenneTwister
    using CompetingClocks: FirstToFire
    using ..VectorAddition
    using ..SampleVAS: sample_sir

    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = FirstToFire{Int64,Float64}()

    starting = zeros(Int, 3 * cnt)
    starting[2:cnt] .= 1
    starting[cnt + 1] = 1  # Start with one infected.
    initial_state = vas_initial(vas, starting)

    simulation = VectorAdditionFSM(vas, initial_state, sampler, rng)
    next_transition = -1
    step_cnt = 0
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
        step_cnt += 1
    end

    # Extract final state
    state = simulation.state.state
    final_S = sum(state[1:cnt])
    final_I = sum(state[cnt+1:2*cnt])
    final_R = sum(state[2*cnt+1:3*cnt])

    @test final_S + final_I + final_R == cnt
    @test final_I == 0
    @test final_R >= 1
    @test step_cnt == (cnt - final_S) + final_R
end
