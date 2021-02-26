using SafeTestsets
using Test

module SampleVAS
function sample_transitions()
    [
        -1  0  1  0  0;
         1 -1 -1 -1  1;
         0  1  0  1 -1
    ]
end
end


@safetestset vas_creates = "VAS creates" begin
using Fleck: VectorAdditionSystem
using ..SampleVAS: sample_transitions
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    @test length(vas.rates) == size(vas.transitions, 2)
end


@safetestset vas_intializes = "VAS initializes state" begin
using Fleck: VectorAdditionSystem, vas_initial, zero_state
using ..SampleVAS: sample_transitions
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    initializer = vas_initial(vas, [1, 1, 0])
    state = zero_state(vas)
    initializer(state)
    @test state == [1, 1, 0]
end


@safetestset vas_reports_hazards = "VAS reports_hazards" begin
using Fleck: VectorAdditionSystem, vas_initial
using Fleck: zero_state, push!, fire!
using ..SampleVAS: sample_transitions
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    initializer = vas_initial(vas, [1, 1, 0])
    state = zero_state(vas)
    disabled = zeros(Int, 0)
    enabled = zeros(Int, 0)
    track_hazards = (idx, dist, enable, rng) -> begin
        if enable == :Enabled
            push!(enabled, idx)
        elseif enable == :Disabled
            push!(disabled, idx)
        else
            @test enable in [:Enabled, :Disabled]
        end
    end
    fire!(track_hazards, vas, state, initializer, "rng")
    @test enabled == [1, 2, 3, 4]
    @test length(disabled) == 0
end



@safetestset vas_changes_state = "VAS hazards during state change" begin
using Fleck: VectorAdditionSystem, vas_input, vas_initial
using Fleck: zero_state, push!, fire!
using ..SampleVAS: sample_transitions
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    initializer = vas_initial(vas, [1, 1, 0])
    state = zero_state(vas)
    initializer(state)
    disabled = zeros(Int, 0)
    enabled = zeros(Int, 0)
    track_hazards = (idx, dist, enable, rng) -> begin
        if enable == :Enabled
            push!(enabled, idx)
        elseif enable == :Disabled
            push!(disabled, idx)
        else
            @test enable in [:Enabled, :Disabled]
        end
    end
    fire_index = 2
    input_change = vas_input(vas, fire_index)
    fire!(track_hazards, vas, state, input_change, "rng")
    @test enabled == [5]
    @test disabled == [2, 3, 4]
end


@safetestset vas_loops = "VAS main loop" begin
using Fleck: VectorAdditionSystem, vas_input, vas_initial
using Fleck: zero_state, push!, fire!
using Random: MersenneTwister
using ..SampleVAS: sample_transitions
    rng = MersenneTwister(2930472)
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    state = zero_state(vas)
    next_transition = nothing
    disabled = zeros(Int, 0)
    enabled = zeros(Int, 0)
    newly_enabled = zeros(Int, 0)
    track_hazards = (idx, dist, enable, rng) -> begin
        if enable == :Enabled
            push!(newly_enabled, idx)
        elseif enable == :Disabled
            push!(disabled, idx)
        else
            @test enable in [:Enabled, :Disabled]
        end
    end
    for i in 1:10
        if isnothing(next_transition)
            input_change = vas_initial(vas, [1, 1, 0])
        else
            input_change = vas_input(vas, next_transition)
        end
        fire!(track_hazards, vas, state, input_change, rng)
        input_change(state)
        if any(state .< 0)
            println("idx $(i) state $(state)")
        end
        @test all(state .>= 0)

        filter!(x -> x ∉ disabled, enabled)
        append!(enabled, newly_enabled)

        if length(enabled) > 0
            next_transition = rand(rng, enabled)
        else
            println("Nothing enabled, $(state)")
            next_transition = nothing
        end
        empty!(newly_enabled)
        empty!(disabled)
    end
end


@safetestset vas_fsm_init = "VAS finite state init" begin
    using Fleck: VectorAdditionFSM, vas_initial, VectorAdditionSystem
    using Fleck: VectorAdditionModel, VectorAdditionFSM, zero_state
    using Fleck: DirectCall, simstep!, vas_input
    using Random: MersenneTwister
    using ..SampleVAS: sample_transitions
    rng = MersenneTwister(2930472)
    transitions = sample_transitions()
    vas = VectorAdditionSystem(transitions, [.5, .8, .7, .7, .7])
    state = zero_state(vas)
    vam = VectorAdditionModel(vas, state, 0.0)
    fsm = VectorAdditionFSM(vam, DirectCall(Int))
    initial_state = vas_initial(vas, [1, 1, 0])
    when, next_transition = simstep!(fsm, initial_state, rng)
    limit = 10
    while next_transition !== nothing && limit > 0
        token = vas_input(vas, next_transition)
        when, next_transition = simstep!(fsm, token, rng)
        limit -= 1
    end
end
