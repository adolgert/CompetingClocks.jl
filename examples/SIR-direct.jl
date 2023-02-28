using Random
using Plots
using Distributions
using Fleck

function sir_vas(β, c, γ)
    take = [
        1 0;
        1 1;
        0 0;
    ]
    give = [
        0 0;
        2 0;
        0 1;
    ]
    rates = [(state) -> Exponential(β*c*state[2]/sum(state)*state[1]), (state) -> Exponential(state[2] * γ)]
    (take, give, rates)
end

p = [0.05,10.0,0.25]; # β,c,γ
u0 = [990,10,0]; # S,I,R

take, give, rates = sir_vas(p...)

rng = MersenneTwister()

vas = VectorAdditionSystem(take, give, rates)
state = zero_state(vas)

vam = VectorAdditionModel(vas, state, 0.0)
fsm = VectorAdditionFSM(vam, DirectCall(Int))

initial_state = vas_initial(vas, u0)

out = Matrix{Float64}(undef, u0[2] + 2*u0[1] + 1, 4)

when, next_transition = simstep!(fsm, initial_state, rng)

event_cnt = 0
tnow = 0.0
out[event_cnt + 1, 2:4] = u0
out[event_cnt + 1, 1] = tnow

while next_transition !== nothing
    tnow += when
    out[event_cnt + 1, 1] = tnow
    out[event_cnt + 1, 2:4] = vam.state
    token = vas_delta(vas, next_transition)
    when, next_transition = simstep!(fsm, token, rng)
    event_cnt += 1
end

plot(
    out[1:event_cnt, 1],
    out[1:event_cnt, 2:4],
    label=["S" "I" "R"],
    xlabel="Time",
    ylabel="Number"
)