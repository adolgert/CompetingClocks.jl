using Random



function sir_vas(N, β, γ)
    take = [
        1 0;
        0 1;
        0 0;
    ]
    give = [
        0 0;
        1 0;
        0 1;
    ]
    rates = [(state) -> state[1] * state[2] * β/N, (state) -> state[2] * γ]
    (take, give, rates)
end

take, give, rates = sir_vas(1000, 0.05, 0.25)

u0 = [990,10,0]; # S,I,R

rng = MersenneTwister()

vas = VectorAdditionSystem(take, give, rates)
state = zero_state(vas)

vam = VectorAdditionModel(vas, state, 0.0)
fsm = VectorAdditionFSM(vam, DirectCall(Int))

initial_state = vas_initial(vas, u0)

out = Matrix{Int64}(undef, u0[2] + 2*u0[1] + 1, 3)

when, next_transition = simstep!(fsm, initial_state, rng)

event_cnt = 0
out[event_cnt + 1, :] = u0

while next_transition !== nothing
    token = vas_delta(vas, next_transition)
    when, next_transition = simstep!(fsm, token, rng)
    event_cnt += 1
    out[event_cnt + 1, :] = vam.state
end
