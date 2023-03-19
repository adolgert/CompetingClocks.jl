using SafeTestsets
using Test


@safetestset vas_sample = "VAS samples with direct" begin
using Random: MersenneTwister
using Fleck: VectorAdditionSystem, VectorAdditionFSM, DirectCall, vas_initial
using Fleck: simstep!, FirstToFire
using ..SampleVAS: sample_sir

    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = DirectCall(Int)

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


@safetestset vas_sample = "VAS samples with first reaction" begin
    using Random: MersenneTwister
    using Fleck: VectorAdditionSystem, VectorAdditionFSM, FirstReaction, vas_initial
    using Fleck: simstep!
    using ..SampleVAS: sample_sir
    
    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = FirstReaction{Int}()

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



@safetestset vas_sample = "VAS samples with first to fire" begin
    using Random: MersenneTwister
    using Fleck: VectorAdditionSystem, VectorAdditionFSM, FirstToFire, vas_initial
    using Fleck: simstep!
    using ..SampleVAS: sample_sir
    
    rng = MersenneTwister(2930472)

    cnt = 30
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    sampler = FirstToFire(Int)

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
