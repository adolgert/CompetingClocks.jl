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
    total_hazard = zero(Float64)
    numerator = zero(Float64)
    for (clock, dist_state) in distributions
        offset = clock == fired_clock ? one(Float64) : zero(Float64)
        fired_hazard = hazard(dist_state.d, when - dist_state.enabling_time)
        total_hazard += fired_hazard
        numerator += (fired_hazard - offset)^2
    end
    return numerator / total_hazard^2
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
