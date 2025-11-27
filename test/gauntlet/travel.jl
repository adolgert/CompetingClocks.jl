module TravelModel

using CompetingClocks
using Random
using Distributions
using Graphs
using EnumX


export travel_commands, travel_make_run, travel_make_graph
export Travel, travel_init_state, travel_enabled
export travel_run, TravelMemory, TravelGraph, TravelRateDist, TravelRateCount
export TravelRateDelay, TravelConfig

@enumx TravelMemory forget remember
@enumx TravelGraph path cycle complete clique
@enumx TravelRateDist exponential general
@enumx TravelRateCount destination pair
@enumx TravelRateDelay none right left some


function travel_make_graph(i, n)
    g = Dict(
            TravelGraph.path => path_graph, # 1
            TravelGraph.cycle => cycle_graph, # 2
            TravelGraph.complete => complete_graph, # 3
            TravelGraph.clique => m-> clique_graph(4, m) # 4
        )[i](n)
    return g
end


# A rate is a `destination` rate when there is one rate per destination location.
function travel_rates_exponential_destination(n, rng)
    hazards = exp.(range(-2, stop=2, length=n))
    return Dict((i, j) => Exponential(inv(hazards[j])) for i in 1:n for j in 1:n if i != j)
end


# A rate is a `travel` rate when there is one rate for each start-finish pair of locations.
function travel_rates_exponential_pair(n, rng)
    hazards = exp.(range(-2, stop=2, length=n*(n-1)))
    rates = Dict{Tuple{Int,Int},Tuple{UnivariateDistribution}}()
    idx = 1
    for i in 1:n, j in 1:n
        i == j && continue
        rates[(i, j)] = Exponential(inv(hazards[idx]))
        idx += 1
    end
    return rates
end


function random_distribution(rng)
    didx = rand(rng, 1:3)
    if didx == 1
        r = -2 + 4 * rand(rng)
        Exponential(inv(exp(r)))
    elseif didx == 2
        alpha = rand(rng, Uniform(0.7, 1.3))
        theta = rand(rng, Uniform(-2, 2))
        Weibull(alpha, inv(exp(theta)))
    elseif didx == 3
        alpha = rand(rng, Uniform(1.0, 6.0))
        theta = rand(rng, Uniform(-2, 2))
        Gamma(alpha, inv(exp(theta)))
    else
        error("Can't make distribution for $didx")
    end
end


function travel_rates_general_destination(n, rng)
    rates = Dict{Tuple{Int,Int},Tuple{UnivariateDistribution}}()
    for target in 1:n
        dist = random_distribution(rng)
        for source in 1:n
            source == target && continue
            rates[(source, target)] = dist
        end
    end
    return rates
end


function travel_rates_general_pair(n, rng)
    rates = Dict{Tuple{Int,Int},Tuple{UnivariateDistribution}}()
    for source in 1:n, target in 1:n
        source == target && continue
        dist = random_distribution(rng)
        rates[(source, target)] = dist
    end
    return rates
end


struct Travel
    g::SimpleGraph{Int64}
    rates::Dict{Tuple{Int,Int},Tuple{UnivariateDistribution,Float64}}
    memory::Vector{Float64}
    remember::TravelMemory.T
end


function draw_delay(delay::TravelRateDelay.T, rng)
    if delay == TravelRateDelay.none
        return zero(Float64)
    elseif delay == TravelRateDelay.left
        return rand(rng, Uniform(-0.5, -0.1))
    elseif delay == TravelRateDelay.right
        return rand(rng, Uniform(0.1, 2.0))
    else
        return rand(rng, Uniform(-0.5, 2.0))
    end
end


function delay_generator(state_cnt::Int, delay::TravelRateDelay.T, count::TravelRateCount.T, rng::AbstractRNG)
    delays = Dict{Tuple{Int,Int},Float64}()
    if count == TravelRateCount.destination
        for destination in 1:state_cnt
            delay_val = draw_delay(delay, rng)
            for source in 1:state_cnt
                source == destination && continue
                delays[(source, destination)] = delay_val
            end
        end
    else
        for destination in 1:state_cnt
            for source in 1:state_cnt
                source == destination && continue
                delay_val = draw_delay(delay, rng)
                delays[(source, destination)] = delay_val
            end
        end
    end
    return delays
end


struct TravelConfig
    memory::TravelMemory.T
    graph::TravelGraph.T
    dist::TravelRateDist.T
    count::TravelRateCount.T
    delay::TravelRateDelay.T
end


function Travel(state_cnt::Int, config::TravelConfig, rng::AbstractRNG)
    g = travel_make_graph(config.graph, state_cnt)
    rates = travel_rates_exponential_destination(state_cnt, rng)
    delays = delay_generator(state_cnt, config.delay, config.count, rng)
    full_rate = Dict((i, j) => (rates[(i, j)], delays[(i, j)]) for (i, j) in keys(rates))
    return Travel(g, full_rate, zeros(Float64, state_cnt), config.memory)
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


function travel_commands(step_cnt::Int, state_cnt::Int, config::TravelConfig, rng::AbstractRNG)
    model = Travel(state_cnt, config, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(step_cnt, sampler, model, rng)
    return commands
end

end
