include("jointreliability.jl")


const nwt_trials = [
    (0.1, true, .25, .26),
    (0.1, false, .25, .26),
    (0.255, true, .255, .26),
    (0.255, false, 1.25, 1.26),
    (0.7, true, 1.25, 1.26),
    (0.7, false, 1.25, 1.26)
]
for (hour, need, ansmin, ansmax) in nwt_trials
    for day in [0.0, 1.0, -5.0, 17.0]
        retmin, retmax = next_work_time(day + hour, 0.25, 0.26, need)
        if abs(retmin - ansmin - day) > 1e-6 || abs(retmax - ansmax - day) > 1e-6
            println("failed on $(hour), $(need) $(retmin) $(retmax)")
        end
    end
end

work_rate = Gamma(9.0, 0.5)
break_rate = LogNormal(3.3, 0.4)
repair_rate = Weibull(1.0, 5.0)
joe = Individual(work_rate, break_rate, repair_rate)

rng = Xoshiro(9378424)
exp = Experiment(20, 10, rng)
