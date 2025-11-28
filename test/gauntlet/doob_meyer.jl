using HypothesisTests


struct DistributionIntegratedHazard
    distributions::Vector{DistributionState}
end

function (dsf::DistributionIntegratedHazard)(t)
    total = zero(Float64)
    for dist_state in dsf.distributions
        total += logccdf(dist_state.d, t - dist_state.enabling_time)
    end
    return exp(total)
end


function doob_meyer(times::Vector{Float64}, distributions::Vector{DistributionState}, when::Float64)
    Gamma = DistributionIntegratedHazard(distributions)
    Gamma_0 = Gamma(when)
    uniform_draws = similar(times)
    for idx in eachindex(times)
        uniform_draws[idx] = 1.0 - exp(-Gamma(times[idx]) + Gamma_0)
    end
    test = ApproximateOneSampleKSTest(uniform_draws, Uniform(0, 1))
    p = pvalue(test)
    (; pvalue=p, supremum_of_difference=test.δ, test=test)
end


function doob_meyer(draws::Vector{ClockDraw}, distributions::Vector{DistributionState}, when::Float64)
    forget_clock = [x[2] for x in draws]
    return doob_meyer(forget_clock, distributions, when)
end
