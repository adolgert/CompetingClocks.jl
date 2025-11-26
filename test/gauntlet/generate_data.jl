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


"""
    parallel_replay(commands, sampler, replica_cnt::Int, rng::Vector)

Returns (vector of samplers, time of last event). The samplers are all initialized
to the last event at the same time with the same previously-fired events.
"""
function parallel_replay(commands, sampler::S, replica_cnt, rng::Vector{T}) where {S <: SSA, T<:AbstractRNG}
    samplers = Vector{S}(undef, replica_cnt)
    for construct_idx in 1:replica_cnt
        samplers[construct_idx] = clone(sampler)
    end
    @threads for run_idx in 1:replica_cnt
        tid = Threads.threadid()
        replay_commands(commands, samplers[run_idx], rng[tid])
    end
    when = commands[end][end]
    return samplers, when
end


"""
Represents a shifted distribution, so it's the distribution and a left-right
shift to the zero-point of the distribution.
"""
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


"""
Represents one sample of one sampler, so it's a clock and a time,
where the clock is an integer in 1:N and the time is greater than the last event
time.
"""
const ClockDraw = Tuple{Int,Float64}


"""
    sample_samplers(samplers, when, rng::Vector)

Samples one event from each sampler. Returns `Vector{ClockDraw}`.
"""
function sample_samplers(samplers, when, rng::Vector{T}) where {T<:AbstractRNG}
    data = similar(samplers, ClockDraw)
    @threads for run_idx in eachindex(samplers)
        tid = Threads.threadid()
        next_when, which = next(samplers[run_idx], when, rng[tid])
        data[run_idx] = (which, next_when)
    end
    return data
end


"""
    retrieve_draws(commands, sampler_cnt, rng)

Given a set of commands, creates `sampler_cnt` samplers and gets one sample
from each of them. Returns `(Vector{ClockDraw}, final_time::Float64)`.
"""
function retrieve_draws(commands, sampler, sampler_cnt, rng)
    samplers, final_time = parallel_replay(commands, sampler, sampler_cnt, rng)
    draws = sample_samplers(samplers, final_time, rng)
    return draws, final_time
end


"""
    jumble!(sample_data::Vector{ClockDraw}, rng)

Used for permutation testing, this mixes which clock goes with which
sampling time.
"""
function jumble!(sample_data::Vector{ClockDraw}, rng)
    clocks = [x[1] for x in sample_data]
    shuffle!(rng, clocks)
    for i in eachindex(clocks)
        sample_data[i] = (clocks[i], sample_data[i][2])
    end
end


function travel_make_run()
    rng = Xoshiro(98327423)
    config = TravelConfig(
        TravelMemory.forget, TravelGraph.path, TravelRateDist.exponential,
        TravelRateCount.destination, TravelRateDelay.none
        )
    model = Travel(2, config, rng)
    sampler = FirstReaction{Int,Float64}()
    commands = travel_run(5, sampler, model, rng)
    sampler2 = FirstReaction{Int,Float64}()
    replay_commands(commands, sampler2, rng)
end
