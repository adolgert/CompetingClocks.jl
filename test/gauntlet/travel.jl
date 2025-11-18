module TravelModel

using CompetingClocks
using Random
using Distributions
using Graphs
using EnumX


export travel_commands, travel_make_run, travel_make_graph, travel_make_rate
export travel_rates_exponential, Travel, travel_init_state, travel_enabled
export travel_run, TravelMemory, TravelGraph

@enumx TravelMemory forget remember
@enumx TravelGraph path cycle complete clique

function travel_make_graph(i, n)
    g = Dict(
            TravelGraph.path => path_graph, # 1
            TravelGraph.cycle => cycle_graph, # 2
            TravelGraph.complete => complete_graph, # 3
            TravelGraph.clique => m-> clique_graph(4, m) # 4
        )[i](n)
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

struct Travel
    g::SimpleGraph{Int64}
    rates::Dict{Tuple{Int,Int},Tuple{UnivariateDistribution,Float64}}
    memory::Vector{Float64}
    remember::TravelMemory.T
end

function Travel(state_cnt::Int, graph_type::TravelGraph.T, memory_type::TravelMemory.T, rng::AbstractRNG)
    g = travel_make_graph(graph_type, state_cnt)
    rates = travel_rates_exponential(state_cnt, rng)
    return Travel(g, rates, zeros(Float64, state_cnt), memory_type)
end
travel_init_state(tr::Travel, rng) = rand(rng, 1:nv(tr.g))
travel_enabled(tr::Travel, state) = Set(neighbors(tr.g, state))

# This sampler is a low-level sampler. This is for testing individual samplers,
# not the SamplingContext. As a result, we track the RNG and the current time.
function travel_run(step_cnt, sampler::SSA, travel_model, rng)
    commands = Vector{Tuple}()
    system_state = travel_init_state(travel_model, rng)
    system_time = 0.0
    enabled_clocks = travel_enabled(travel_model, system_state)
    for clock in enabled_clocks
        dist, offset = travel_model.rates[(system_state, clock)]
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
            if travel_model.remember == TravelMemory.remember
                offset -= travel_model.memory[en_idx]
            end
            enable!(sampler, en_idx, en_rate, when + offset, when, rng)
            push!(commands, (:enable, en_idx, en_rate, when + offset, when))
        end
        for un_idx in intersect(enabled_clocks, next_enabled)
            # Even if a clock was enabled before, the model might define a rate from
            # state 2->3 that's different from the rate from state 7->3. In that case,
            # disable and enable.
            previous_rate, previous_offset = travel_model.rates[(system_state, un_idx)]
            current_rate, current_offset = travel_model.rates[(which, un_idx)]
            if previous_rate != current_rate || previous_offset != current_offset
                # Flip a coin on whether we silently re-enable. Either way should
                # be fine for all samplers, so let's check that it's fine.
                if rand(rng, Bool)
                    disable!(sampler, un_idx, when)
                    push!(commands, (:disable, un_idx, when))
                    enable!(sampler, un_idx, current_rate, when + current_offset, when, rng)
                    push!(commands, (:enable, un_idx, current_rate, when + current_offset, when))
                else
                    enable!(sampler, un_idx, current_rate, when + current_offset, when, rng)
                    push!(commands, (:enable, un_idx, current_rate, when + current_offset, when))
                end
            end
        end
        enabled_clocks = next_enabled
        system_time = when
        system_state = which
        when, which = next(sampler, system_time, rng)
    end
    return commands
end


function travel_commands(step_cnt, state_cnt::Int, graph_type::TravelGraph.T, memory_type::TravelMemory.T, rng::AbstractRNG)
    model = Travel(state_cnt, graph_type, memory_type, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(step_cnt, sampler, model, rng)
    return commands
end

end
