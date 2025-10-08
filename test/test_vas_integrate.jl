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
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
    end
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
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
    end
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
    while next_transition !== nothing
        when, next_transition = simstep!(simulation)
    end
end
