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
    rates = [(state) -> Exponential(1.0/(β*c*state[2]/sum(state)*state[1])), (state) -> Exponential(1.0/(state[2] * γ))]
    (take, give, rates)
end

p = [0.05,10.0,0.25]; # β,c,γ
u0 = [990,10,0]; # S,I,R

take, give, rates = sir_vas(p...)

rng = MersenneTwister()
vas = VectorAdditionSystem(take, give, rates)

# sampler = DirectCall{Int}()
# sampler = FirstReaction{Int}()
sampler = NextReaction{Int}()
# sampler = FirstToFire{Int}()
fsm = VectorAdditionFSM(vas, vas_initial(vas, u0), sampler, rng)
out = Matrix{Float64}(undef, u0[2] + 2*u0[1] + 1, 4)

event_cnt = 0
while true
    when, next_transition = simstep!(fsm)
    if next_transition === nothing
        break
    end
    event_cnt += 1
    out[event_cnt, 1] = fsm.state.when
    out[event_cnt, 2:4] = fsm.state.state
end

plot(
    out[1:event_cnt, 1],
    out[1:event_cnt, 2:4],
    label=["S" "I" "R"],
    xlabel="Time",
    ylabel="Number"
)