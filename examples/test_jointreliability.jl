include("jointreliability.jl")

using Logging
using ColorSchemes
using StatsPlots

for hour in [0.005, 0.7]
    for day in [0.0, 1.0, -5.0, 17.0]
        now = day + hour
        retmin, retmax = next_work_time(now, 0.01)
        @assert retmin >= 0.0
        @assert retmin < 0.99
        @assert retmax > retmin
        @assert retmax - retmin <= 0.01 + 1e-6
    end
end

work_rate = LogUniform(.5, 0.99) # Gamma(9.0, 0.5)
break_rate = LogNormal(3.3, 0.4)
repair_rate = Weibull(1.0, 5.0)
joe = Individual(work_rate, break_rate, repair_rate)

struct OnlyMyLogger <: AbstractLogger
    instance::ConsoleLogger
    OnlyMyLogger() = new(ConsoleLogger(stdout, Logging.Debug))
end
function Logging.shouldlog(logger::OnlyMyLogger, level, _module, group, id)
    _module == Main && Logging.shouldlog(logger.instance, level, _module, group, id)
end

Logging.min_enabled_level(logger::OnlyMyLogger) = Logging.min_enabled_level(logger.instance)
Logging.handle_message(logger::OnlyMyLogger, args...; kwargs...) = Logging.handle_message(logger.instance, args...; kwargs...)

function plotgrid(observation::ObserveHistogram)
    heatmap(observation.counts, c = cgrad(:acton))
    xlabel!("Broken")
    ylabel!("Worked")
    png("counts.png")
end    

function compare_across_workers(obs::Vector{ObserveHistogram}, labels, title)
    firstplot = true
    for (obs_idx, observation) in enumerate(obs)
        worker_cnt = total_workers(observation)
        broken = vec(sum(observation.counts, dims=1))
        cnt = findlast(broken .> 0)
        normed = broken[1:cnt] / sum(broken)
        if firstplot
            firstplot = false
            plot(1:cnt, normed, seriestype=:scatter, label=labels[obs_idx], markersize=8)
        else
            plot!(1:cnt, normed, seriestype=:scatter, label=labels[obs_idx], markersize=8)
        end
    end
    xlabel!("Count of Broken")
    ylabel!("Probability Mass Function")
    title!(title)
    png(replace(title, " "=>"") * ".png")
end


function plot_timeline(obs::ObserveContinuous, experiment::Experiment)
    state = obs.state
    last_time = state[end].time
    first_idx = findlast([x.time < last_time - 5 for x in state])
    times = [x.time for x in state[first_idx:end]]
    times .-= 4009
    working = [x.working for x in state[first_idx:end]]
    broken = [x.broken for x in state[first_idx:end]]
    ready = [worker_cnt(experiment) - x.working - x.broken for x in state[first_idx:end]]
    plot(times, working, label="working", line=(:steppost, 2))
    plot!(times, broken, label="broken", line=(:steppost, 2))
    # plot!(times, ready, label="ready", line=(:steppost, 2))
    xlabel!("Time [days]")
    ylabel!("Status [count]")
    title!("Timeline of Work Crew Over Five Days")
    png("timeline_month.png")
end


function show_number_broken()
    rng = Xoshiro(9378424)
    observations = ObserveHistogram[]
    labels = String[]
    years = 100
    day_cnt = 365 * years

    for worker_cnt in [10, 11, 12, 13, 14, 16, 20, 100]
        experiment = Experiment(worker_cnt, 10, rng)
        burn = min(day_cnt ÷ 10, 3650)
        observation = ObserveHistogram(experiment, burn)
        run(experiment, observation, day_cnt)
        push!(observations, observation)
        push!(labels, string(worker_cnt))
    end

    compare_across_workers(observations, labels, "Number Broken as Workers Increase")
end


function walk_simulation()
    debug_logger = OnlyMyLogger()
    rng = Xoshiro(234987)
    with_logger(debug_logger) do
        day_cnt = 20
        experiment = Experiment(worker_cnt, 10, rng)
        observation = ObserveLots(day_cnt, worker_cnt(experiment))
        run(experiment, observation, day_cnt)
        plot(1:day_cnt, observation.status[2, :], seriestype=:scatter, label="repair")
        plot!(1:day_cnt, observation.status[1, :], seriestype=:scatter, label="working")
        png("record.png")
        println(100 * observation.broken_duration / day_cnt)
    end
end


function show_competition_effect()
    rng = Xoshiro(4377124)
    observations = ObserveHistogram[]
    labels = String[]
    years = 1000
    day_cnt = 365 * years
    worker_cnt = 20
    for must_work in [1, 5, 10, 15, 20]
        experiment = Experiment(worker_cnt, must_work, rng)
        burn = min(day_cnt ÷ 10, 3650)
        observation = ObserveHistogram(experiment, burn)
        run(experiment, observation, day_cnt)
        push!(observations, observation)
        push!(labels, string(must_work))
    end

    compare_across_workers(observations, labels, "Number Broken as Crew Increases")
end


function show_typical_timeline()
    rng = Xoshiro(9234232)
    years = 11
    day_cnt = 365 * years
    worker_cnt = 16
    observation = ObserveContinuous
    experiment = Experiment(worker_cnt, 10, rng)
    observation = ObserveContinuous()
    run(experiment, observation, day_cnt)

    plot_timeline(observation, experiment)
end

# show_competition_effect()
# show_number_broken()
# show_typical_timeline()
