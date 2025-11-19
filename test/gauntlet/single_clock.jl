using Base.Threads: @threads, maxthreadid, threadid
using CompetingClocks: hazard

"""
    mark_calibration_brier(distributions, fired_clock, when)

Calculating one component of the Brier score. The whole Brier score is

```math
    \\frac{1}{N}\\sum_i^N \\sum_i^R (f_{ti}-o_{ti})^2
```

where ``o_{ti}`` is 1 for the fired clock and 0 otherwise. Here ``f_{ti}`` is
the hazard of one clock divided by the total hazard. This function calculates
the above for one outcome ``i``.
"""
function mark_calibration_brier(distributions, fired_clock, when)
    total_hazard = 0.0
    hazards = Dict{Int,Float64}()
    for (clock, ds) in distributions
        h = hazard(ds.d, when - ds.enabling_time)
        hazards[clock] = h
        total_hazard += h
    end
    total_hazard == 0.0 && return 0.0

    s = 0.0
    for (clock, h) in hazards
        p = h / total_hazard
        y = clock == fired_clock ? 1.0 : 0.0
        s += (p - y)^2
    end
    # Divide by distributions because this it the number of classes
    # for multiclass Brier.
    return s / length(distributions)
end



"""
    mark_calibration_conditional_time(draws, distributions)

This tests whether ``P[E|T]`` is correct. It's the probability that the right
clock fires given when a clock fired. It's a simple value because it's the hazard
for the clock that fired divided by the sum of all hazards at the time that clock
fired.
"""
function mark_calibration_conditional_time(draws, distributions)
    # Use thread-local storage for reduction
    partial_sums = zeros(Float64, Threads.maxthreadid())

    @threads for i in eachindex(draws)
        tid = Threads.threadid()
        partial_sums[tid] += mark_calibration_brier(distributions, draws[i][1], draws[i][2])
    end

    return sum(partial_sums) / length(draws)
end
