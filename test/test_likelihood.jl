using SafeTestsets

module LikelihoodHelper
using CompetingClocks
using Distributions
using Test
export Fire, Enable, Disable, LikelihoodTestCase, ClockSpec, ClockState, execute
export SampleState
abstract type LHAction end
struct Fire <: LHAction
    clock::Symbol
    time::Float64
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
    likelihood::Vector{Float64}
    ClockState() = new(Dict{Symbol,ClockSpec}(), 0.0, Float64[])
end
function execute(cs::ClockState, action::Fire)
    @test action.clock ∈ keys(cs.clocks)
    @test cs.clocks[action.clock].enabled
    push!(cs.likelihood, step_likelihood(cs, action.clock, action.time))
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
struct SampleState
    sampler::Any
    step_likelihood::Vector{Float64}
    path_likelihood::Vector{Float64}
    SampleState(sampler) = new(sampler, Float64[], Float64[])
end
function execute(cs::SampleState, action::Fire)
    push!(cs.step_likelihood, steploglikelihood(cs.sampler, action.clock))
    fire!(cs.sampler, action.clock, action.time)
    push!(cs.path_likelihood, trajectoryloglikelihood(cs.sampler, time(cs.sampler)))
    cs
end
execute(cs::SampleState, action::Enable) = (enable!(cs.sampler, action.clock, action.distribution, action.offset); cs)
execute(cs::SampleState, action::Disable) = (disable!(cs.sampler, action.clock); cs)

function execute(cs::ClockState, sampler::SampleState, action::Vector{LHAction})
    for idx in eachindex(action)
        execute(cs, action[idx])
        execute(sampler, action[idx])
        if action[idx] isa Fire
            @test length(cs.likelihood) == length(sampler.step_likelihood)
            @test length(cs.likelihood) == length(sampler.path_likelihood)
            @test abs(cs.likelihood[end] - sampler.step_likelihood[end]) < 1e-6
            @test abs(sum(cs.likelihood) - sampler.path_likelihood[end]) < 1e-6
        end
    end
end

function step_likelihood(cs::ClockState, fired_idx, t1)
    # t0 was the last firing time. t1 is the next firing time.
    # fired_idx fires at time t1.
    t0 = cs.time # Calculating the likelihood from t0 to t1.
    # spec.t0 is simulation time when distribution was enabled.
    # spec.te is absolute time for the zero-point of the distribution.
    ll = 0.0
    @test t0 < t1  # Assert invariant that time steps aren't simultaneous.
    for spec in values(cs.clocks)
        @test spec.t0 <= t0  # Clock can't be enabled after last time point.
        !spec.enabled && continue
        te = spec.te
        if spec.t0 == te  # enabled at the zero-point of the distribution.
            if spec.clock == fired_idx
                @test te < t1  # Can't fire unless pdf is nonzero.
                ll += logpdf(spec.distribution, t1 - te)
            else
                if t1 > spec.te
                    ll += logccdf(spec.distribution, t1 - te)
                end
            end
            if t0 > te
                ll -= logccdf(spec.distribution, t0 - te)
            end
        elseif te < spec.t0 # left-shifted
            if spec.clock == fired_idx
                ll += logpdf(spec.distribution, t1 - te)
            else
                ll += logccdf(spec.distribution, t1 - te)
            end
            ll -= logccdf(spec.distribution, t0 - te)
        elseif spec.t0 < te && te < t1 # right-shifted
            if spec.clock == fired_idx
                @test te < t1  # Can't fire unless pdf is nonzero.
                ll += logpdf(spec.distribution, t1 - te)
            else
                if t1 > te
                    ll += logccdf(spec.distribution, t1 - te)
                end
            end
        elseif spec.t0 < te && te >= t1 # right-shifted
            nothing  # No contribution in this case. Wasn't fireable.
        else
            error("Can't get here.")
        end
    end
    return ll
end

end # module LikelihoodHelper

@safetestset likely_single_clock = "Single clock" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end
