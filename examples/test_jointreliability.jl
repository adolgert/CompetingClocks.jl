include("jointreliability.jl")

using Logging

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

rng = Xoshiro(9378424)
experiment = Experiment(16, 10, rng)

plot(experiment.group[1].work_dist)
title!("Work Duration Distribution")
png("work.png")
plot(experiment.group[1].fail_dist)
title!("Distribution of Failures")
png("fail.png")
plot(experiment.group[1].repair_dist)
title!("Repair Time Durations")
png("fix.png")


struct OnlyMyLogger <: AbstractLogger
    instance::ConsoleLogger
    OnlyMyLogger() = new(ConsoleLogger(stdout, Logging.Debug))
end
function Logging.shouldlog(logger::OnlyMyLogger, level, _module, group, id)
    _module == Main && Logging.shouldlog(logger.instance, level, _module, group, id)
end

Logging.min_enabled_level(logger::OnlyMyLogger) = Logging.min_enabled_level(logger.instance)
Logging.handle_message(logger::OnlyMyLogger, args...; kwargs...) = Logging.handle_message(logger.instance, args...; kwargs...)

debug_logger = OnlyMyLogger()
with_logger(debug_logger) do
    day_cnt = 1
    observation = run(experiment, day_cnt)
    day_cnt = days(observation)
    plot(1:day_cnt, observation.status[2, :], seriestype=:scatter, label="repair")
    plot!(1:day_cnt, observation.status[1, :], seriestype=:scatter, label="working")
    png("record.png")
    println(100 * observation.broken_duration / day_cnt)
end
