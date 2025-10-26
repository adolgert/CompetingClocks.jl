module LikelihoodHelper
using Distributions

abstract type LHAction end
struct Fire <: LHAction
    time::Float64
    clock::Symbol
end
struct Enable <: LHAction
    clock::Symbol
    distribution::UnivariateDistribution
    offset::Float64
end
struct Disable <: LHAction
    clock::Symbol
end
struct LikelihoodTestCase
    name::String
    steps::Vector{LHAction}
end
mutable struct ClockSpec
    clock::Symbol
    distribution::UnivariateDistribution
    te::Float64        # Distribution zero-point (absolute time)
    t0::Float64      # When it was enabled (absolute time)
    enabled::Bool
end
mutable struct ClockState
    clocks::Dict{Symbol,ClockSpec}
    time::Float64
end
function execute(cs::ClockState, action::Fire)
    @assert action.clock ∈ keys(cs.clocks)
    @assert cs.clocks[action.clock].enabled
    cs.time = action.time
    cs.clocks[action.clock].enabled = false
    cs
end
function execute(cs::ClockState, action::Enable)
    if action.clock ∈ keys(cs.clocks)
        cs.clocks[action.clock].distribution = action.distribution
        cs.clocks[action.clock].te = cs.time + action.offset
        cs.clocks[action.clock].t0 = cs.time
        cs.clocks[action.clock].enabled = true
    else
        cs.clocks[action.clock] = ClockSpec(
            action.clock, action.distribution, cs.time + action.offset,
            cs.time, true
        )
    end
    cs
end
function execute(cs::ClockState, action::Disable)
    @assert action.clock ∈ keys(cs.clocks)
    cs.clocks[action.clock].enabled = false
    cs
end
function execute(cs::ClockState, action::Vector{LHAction})
    foldl(execute, cs, action)
end
function oracle_step_likelihood(distributions, t0, t1, fired_idx)
    # distributions is a vector of (dist, te, when_enabled) tuples
    ll = 0.0
    for (i, (dist, te, when_enabled)) in enumerate(distributions)
        # Contribution from this distribution
        if i == fired_idx
            ll += logpdf(dist, t1 - te)
        else
            ll += logccdf(dist, t1 - te)
        end
        # Subtract the baseline at t0
        if t0 > te
            ll -= logccdf(dist, t0 - te)
        end
        # Adjust for late enabling
        if when_enabled > te && when_enabled > t0
            ll -= logccdf(dist, when_enabled - te)
        end
    end
    return ll
end


end # module LikelihoodHelper
