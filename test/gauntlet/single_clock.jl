
@raw"""
    mark_calibration_brier(distributions, fired_clock, when)

Calculating one component of the Brier score. The whole Brier score is

```math
    \frac{1}{N}\sum_i^N \sum_i^R (f_{ti}-o_{ti})^2
```

where ``o_{ti}`` is 1 for the fired clock and 0 otherwise. Here ``f_{ti}`` is
the hazard of one clock divided by the total hazard. This function calculates
the above for one outcome ``i``.
"""
function mark_calibration_brier(distributions, fired_clock, when)
    total_hazard = zero(Float64)
    numerator = zero(Float64)
    for (clock, (dist, enabling)) in distributions
        if clock == fired_clock
            fired_hazard = hazard(dist, when - enabling)
            total_hazard += fired_hazard
            numerator += (fired_hazard - 1)^2
        else
            unfired_hazard = hazard(dist, when - enabling)
            total_hazard += unfired_hazard
            numerator += unfired_hazard^2
        end
    end
    return numerator / total_hazard^2
end


function mark_calibration_conditional_time(draws, distributions)
    total = @reduce(+, @threads :static for i in eachindex(draws)
        mark_calibration_brier(distributions, draws[i][1], draws[i][2])
    end)
    return total
end
