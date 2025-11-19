using Base.Threads: @threads

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


function parallel_replay(commands, replica_cnt, rng::Vector{T}) where {T<:AbstractRNG}
    samplers = Vector{FirstReaction{Int,Float64}}(undef, replica_cnt)
    for construct_idx in 1:replica_cnt
        samplers[construct_idx] = FirstReaction{Int,Float64}()
    end
    @threads for run_idx in 1:replica_cnt
        tid = Threads.threadid()
        replay_commands(commands, samplers[run_idx], rng[tid])
    end
    when = commands[end][end]
    return samplers, when
end


struct DistributionState
    d::UnivariateDistribution
    enabling_time::Float64
end


"""
    final_enabled_distributions(commands, dist_cnt)

Returns a dict with enabled distributions at the final event.
Returns a `Dict{Int,DistributionState}()`.
"""
function final_enabled_distributions(commands)
    dist = Dict{Int,Union{DistributionState,Nothing}}()
    for idx in reverse(eachindex(commands))
        cmd = commands[idx]
        if cmd[1] == :disable || cmd[1] == :fire
            dist_idx = cmd[2]
            if !haskey(dist, dist_idx)
                dist[dist_idx] = nothing
            end
        elseif cmd[1] == :enable
            dist_idx = cmd[2]
            if !haskey(dist, dist_idx)
                dist[dist_idx] = DistributionState(cmd[3], cmd[4])
            end
        else
            error("Unknown command $(cmd)")
        end
    end
    return Dict(k => v for (k, v) in dist if !isnothing(v))
end


function sample_samplers(samplers, when, rng::Vector{T}) where {T<:AbstractRNG}
    data = similar(samplers, Tuple{Int,Float64})
    @threads for run_idx in eachindex(samplers)
        tid = Threads.threadid()
        next_when, which = next(samplers[run_idx], when, rng[tid])
        data[run_idx] = (which, next_when)
    end
    return data
end


function jumble!(sample_data::Vector{Tuple{Int,Float64}}, rng)
    reorder = shuffle(rng, 1:length(sample_data))
    clocks = [x[1] for x in sample_data]
    shuffle!(rng, clocks)
    for i in eachindex(clocks)
        sample_data[i] = (clocks[i], sample_data[i][2])
    end
end


function travel_make_run()
    rng = Xoshiro(98327423)
    model = Travel(2, TravelGraph.path, TravelMemory.forget, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng)
    sampler2 = FirstReaction{Int,Float64}()
    replay_commands(commands, sampler2, rng)
end
