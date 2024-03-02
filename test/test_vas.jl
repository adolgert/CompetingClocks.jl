using SafeTestsets
using Test



@safetestset vas_creates = "VAS creates" begin
using ..VectorAddition
using ..SampleVAS: sample_transitions
    take, give, rates = sample_transitions()
    vas = VectorAdditionSystem(take, give, rates)
    @test length(vas.rates) == size(vas.take, 2)
end


@safetestset vas_intializes = "VAS initializes state" begin
using ..VectorAddition
using ..SampleVAS: sample_transitions
    take, give, rates = sample_transitions()
    vas = VectorAdditionSystem(take, give, rates)
    initializer = vas_initial(vas, [1, 1, 0])
    state = zero_state(vas)
    initializer(state)
    @test state == [1, 1, 0]
end


@safetestset vas_reports_hazards = "VAS reports_hazards" begin
using Fleck: push!, DebugWatcher, enable!, disable!
using ..VectorAddition
using ..SampleVAS: sample_transitions
take, give, rates = sample_transitions()
vas = VectorAdditionSystem(take, give, rates)
initializer = vas_initial(vas, [1, 1, 0])

state = zero_state(vas)
track_hazards = DebugWatcher{Int}()
fire!(track_hazards, vas, state, initializer, 0.0, "rng")
enabled = Set(entry.clock for entry in track_hazards.enabled)
@test enabled == Set([1, 2, 3, 4])
@test length(track_hazards.disabled) == 0
end



@safetestset vas_changes_state = "VAS hazards during state change" begin
using ..VectorAddition
using ..SampleVAS: sample_transitions, DebugWatcher, enable!, disable!
take, give, rates = sample_transitions()
vas = VectorAdditionSystem(take, give, rates)
initializer = vas_initial(vas, [1, 1, 0])

state = zero_state(vas)
initializer(state)
track_hazards = DebugWatcher{Int}()
fire_index = 2
input_change = vas_delta(vas, fire_index)
fire!(track_hazards, vas, state, input_change, 0.0, "rng")
enabled = Set(entry.clock for entry in track_hazards.enabled)
@test enabled == Set([5])
disabled = Set(entry.clock for entry in track_hazards.disabled)
@test disabled == Set([2, 3, 4])
end


@safetestset vas_loops = "VAS main loop" begin
# This test takes apart the main loop in order to examine whether
# the model works without a sampler. It watches VAS work.
using Fleck: TrackWatcher, enable!, disable!
using Random: MersenneTwister
using ..VectorAddition
using ..SampleVAS: sample_transitions
    rng = MersenneTwister(2930472)
    take, give, rates = sample_transitions()
    vas = VectorAdditionSystem(take, give, rates)
    state = zero_state(vas)
    next_transition = nothing
    disabled = zeros(Int, 0)
    enabled = zeros(Int, 0)
    newly_enabled = zeros(Int, 0)
    track_hazards = TrackWatcher{Int}()
    curtime = 0.0
    for i in 1:10
        if isnothing(next_transition)
            input_change = vas_initial(vas, [1, 1, 0])
        else
            input_change = vas_delta(vas, next_transition)
        end
        fire!(track_hazards, vas, state, input_change, curtime, rng)
        if any(state .< 0)
            println("idx $(i) state $(state)")
        end
        @test all(state .>= 0)

        if length(track_hazards) > 0
            rand_idx = rand(rng, 1:length(track_hazards))
            select_idx = 1
            for entry in track_hazards
                if select_idx == rand_idx
                    next_transition = entry.clock
                end
                select_idx+=1
            end
            curtime += 1.0
        else
            println("Nothing enabled, $(state)")
            next_transition = nothing
        end
    end
end


@safetestset vas_fsm_init = "VAS finite state init" begin
    using Fleck: DirectCall
    using Random: MersenneTwister
    using ..VectorAddition
    using ..SampleVAS: sample_transitions
    rng = MersenneTwister(2930472)
    take, give, rates = sample_transitions()
    vas = VectorAdditionSystem(take, give, rates)
    initial_state = vas_initial(vas, [1, 1, 0])
    fsm = VectorAdditionFSM(vas, initial_state, DirectCall{Int}(), rng)
    when, next_transition = simstep!(fsm)
    limit = 10
    while next_transition !== nothing && limit > 0
        when, next_transition = simstep!(fsm)
        limit -= 1
    end
end


@safetestset vas_fsm_sir = "VAS runs SIR to completion" begin
    using Fleck: DirectCall
    using Random: MersenneTwister
    using ..VectorAddition
    using ..SampleVAS: sample_sir
    rng = MersenneTwister(979797)
    cnt = 10
    vas = VectorAdditionSystem(sample_sir(cnt)...)
    starting = zeros(Int, 3 * cnt)
    starting[2:cnt] .= 1
    starting[cnt + 1] = 1  # Start with one infected.
    initial_state = vas_initial(vas, starting)
    fsm = VectorAdditionFSM(vas, initial_state, DirectCall{Int}(), rng)
    when, next_transition = simstep!(fsm)
    event_cnt = 0
    while next_transition !== nothing
        when, next_transition = simstep!(fsm)
        event_cnt += 1
    end
    println(event_cnt)
    @test event_cnt > 0
    @test sum(fsm.state.state) == cnt
end
