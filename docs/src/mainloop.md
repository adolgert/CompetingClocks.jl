# Sample Main Loop

Let's walk through a small simulation so that we can see how and when Fleck could appear in the main loop.

Start by defining a state for the simulation.

```julia
mutable struct PhysicalState
    population::Vector{Int}
end

function initialize!(physical::PhysicalState)
    physical.population .= 0
end
```

We need to organize the simulation, itself, which carries the state of the clocks, as well as the physical state. We'll call it a finite state machine (FSM) because it has the traits of a [Moore machine](https://en.wikipedia.org/wiki/Moore_machine).

```julia
mutable struct SimulationFSM
    physical::PhysicalState
    sampler::SSA
    rng::Xoshiro
    initializer::Function
    is_initialized::Bool
    observer::Function
end
```


```julia
function simstep!(state::SimState)
    if !fsm.is_initialized
        fire!(fsm.sampler, fsm.vas, fsm.state.state, fsm.initializer, fsm.state.when, fsm.rng)
        fsm.is_initialized = true
    end
    (when, what) = next(fsm.sampler, fsm.state.when, fsm.rng)
    if when < Inf
        action = vas_delta(fsm.vas, what)
        fsm.state.when = when
        fire!(fsm.sampler, fsm.vas, fsm.state.state, action, fsm.state.when, fsm.rng)
        fsm.observer(fsm.state, when, what)
    else
        (when, what)
    end
end
```


```julia
function fire!(sampler, vas::VectorAdditionSystem, state, modify_state, now, rng)
    former = copy(state)
    modify_state(state)
    for rate_idx in eachindex(vas.rates)
        was_enabled = all(former .- vas.take[:, rate_idx] .>= 0)
        now_enabled =  all(state .- vas.take[:, rate_idx] .>= 0)
        if was_enabled && !now_enabled
            disable!(sampler, rate_idx, now)
        elseif !was_enabled && now_enabled
            enable!(sampler, rate_idx, vas.rates[rate_idx](state), now, now, rng)
        elseif was_enabled && now_enabled
            ratefunc = vas.rates[rate_idx]
            former_rate = ratefunc(former)
            current_rate = ratefunc(state)
            if former_rate != current_rate
                enable!(sampler, rate_idx, current_rate, now, now, rng)
            end  # Else don't notify because rate is the same.
        end
    end
end
```
