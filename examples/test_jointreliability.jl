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

work_rate = Gamma(9.0, 0.5)
break_rate = LogNormal(3.3, 0.4)
repair_rate = Weibull(1.0, 5.0)
joe = Individual(work_rate, break_rate, repair_rate)

rng = Xoshiro(9378424)
experiment = Experiment(20, 10, rng)

plot(experiment.group[1].work_dist)
title!("Work Duration Distribution")
png("work.png")
plot(experiment.group[1].fail_dist)
title!("Distribution of Failures")
png("fail.png")
plot(experiment.group[1].repair_dist)
title!("Repair Time Durations")
png("fix.png")


# debug_logger = ConsoleLogger(stderr, Logging.Debug)
# with_logger(debug_logger) do
    result = run(experiment, 1000)
    println(result)
# end
