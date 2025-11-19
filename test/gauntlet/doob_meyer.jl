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
        uniform_draws[idx] = 1.0 - exp(-Gamma(times[idx]) - Gamma_0)
    end
    ApproximateOneSampleKSTest(uniform_draws, Uniform(0, 1))
end
