using CompetingClocks
using Random
using Distributions
using Graphs


function travel_make_graph(i, n)
    g = [
        path_graph, # 1
        cycle_graph, # 2
        complete_graph, # 3
        m-> clique_graph(4, m) # 4
        ][i](n)
    return g
end


function travel_make_rate(i, rng)
    if i == 1
        # Second member of tuple is a time offset.
        (Exponential(0.2 * rand(rng, Float64) + 1.0), 0.0)
    else
        error("make_rate $i isn't defined")
    end
end

function travel_rates_exponential(n, rng)
    hazards = [travel_make_rate(1, rng) for _ in 1:n]
    return Dict((i, j) => hazards[j] for i in 1:n for j in 1:n)
end

struct TravelModel
    g::SimpleGraph{Int64}
    rates::Dict{Tuple{Int,Int},Tuple{UnivariateDistribution,Float64}}
    memory::Vector{Float64}
    remember::Bool
end
function TravelModel(state_cnt, graph_type, memory_type, rng)
    g = travel_make_graph(graph_type, state_cnt)
    rates = travel_rates_exponential(state_cnt, rng)
    TravelModel(g, rates, zeros(Float64, state_cnt), memory_type == :remember)
end
travel_init_state(tr::TravelModel, rng) = rand(rng, nv(tr.g))
travel_enabled(tr::TravelModel, state) = Set(neighbors(tr.g, state))

# This sampler is a low-level sampler. This is for testing individual samplers,
# not the SamplingContext. As a result, we track the RNG and the current time.
function travel_run(step_cnt, sampler::SSA, travel_model, rng)
    commands = Vector{Tuple}()
    init_state = travel_init_state(travel_model, rng)
    system_time = 0.0
    enabled_clocks = travel_enabled(travel_model, init_state)
    for clock in enabled_clocks
        dist, offset = travel_model.rates[(init_state, clock)]
        enable!(sampler, clock, dist, offset, system_time, rng)
        push!(commands, (:enable, clock, dist, offset, system_time))
    end
    when, which = next(sampler, system_time, rng)
    for _ in 1:step_cnt
        @assert isfinite(when)
        @assert when > system_time
        for mem_idx in enabled_clocks
            travel_model.memory[mem_idx] += when - system_time
        end
        travel_model.memory[which] = 0.0
        fire!(sampler, which, when)
        push!(commands, (:fire, which, when))
        delete!(enabled_clocks, which)
        next_enabled = travel_enabled(travel_model, which)
        to_disable = setdiff(enabled_clocks, next_enabled)
        to_enable = setdiff(next_enabled, enabled_clocks)
        for dis_idx in to_disable
            disable!(sampler, dis_idx, when)
            push!(commands, (:disable, dis_idx, when))
        end
        for en_idx in to_enable
            en_rate, offset = travel_model.rates[(which, en_idx)]
            if travel_model.remember
                offset -= travel_model.memory[en_idx]
            end
            enable!(sampler, en_idx, en_rate, when + offset, when, rng)
            push!(commands, (:enable, en_idx, en_rate, when + offset, when))
        end
        enabled_clocks = next_enabled
        system_time = when
        when, which = next(sampler, system_time, rng)
    end
    return commands
end


function replay_commands(commands, sampler, rng)
    for command in commands
        if command[1] == :enable
            enable!(sampler, command[2:end]..., rng)
        elseif command[1] == :fire
            fire!(sampler, command[2:end]...)
        elseif command[1] == :disable
            disable!(sampler, command[2:end]...)
        else
            error("unknown command $(command[1])")
        end
    end
end

function travel_make_run()
    rng = Xoshiro(98327423)
    model = TravelModel(2, 1, :forget, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng)
    sampler2 = FirstReaction{Int,Float64}()
    replay_commands(commands, sampler2, rng)
end
