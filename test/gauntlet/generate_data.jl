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


function parallel_replay(commands, replica_cnt, base_seed)
    samplers = Vector{FirstReaction{Int,Float64}}(undef, replica_cnt)
    rng = [Xoshiro(base_seed + i) for i in 1:replica_cnt]
    for construct_idx in 1:replica_cnt
        samplers[construct_idx] = FirstReaction{Int,Float64}()
    end
    @threads for run_idx in 1:replica_cnt
        replay_commands(commands, samplers[run_idx], rng[run_idx])
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
function final_enabled_distributions(commands, dist_cnt)
    dist = Vector{Union{DistributionState,Nothing}}(undef, dist_cnt)
    fill!(dist, nothing)
    found = fill(false, dist_cnt)
    for idx in reverse(eachindex(commands))
        cmd = commands[idx]
        if cmd[1] == :disable || cmd[1] == :fire
            dist_idx = cmd[2]
            if !found[dist_idx]
                dist[dist_idx] = nothing
                found[dist_idx] = true
            end
        elseif cmd[1] == :enable
            dist_idx = cmd[2]
            if !found[dist_idx]
                dist[dist_idx] = DistributionState(cmd[3], cmd[4])
                found[dist_idx] = true
            end
        end
    end
    enabled = Dict{Int,DistributionState}()
    for copy_idx in eachindex(dist)
        if !isnothing(dist[copy_idx])
            enabled[copy_idx] = dist[copy_idx]
        end
    end
    return enabled
end


function generate_data(samplers, when, rng)
    data = similar(samplers, Tuple{Int,Float64})
    @threads for run_idx in eachindex(samplers)
        next_when, which = next(samplers[run_idx], when, rng)
        data[run_idx] = (which, next_when)
    end
    return data
end
