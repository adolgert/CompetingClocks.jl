using SafeTestsets

module LikelihoodHelper
using CompetingClocks
using Distributions
using Test
using QuadGK
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
    check_sum::Bool
    ClockState() = new(Dict{Symbol,ClockSpec}(), 0.0, Float64[], false)
end
function core_matrix_sum(cs::ClockState)
    # Integrate the likelihood of a clock to get the marginal probability that clock fires.
    # Sum over all clocks to get normalization, which should be 1.
    marginal = Dict{Symbol,Float64}()
    for clock in keys(cs.clocks)
        if cs.clocks[clock].enabled
            marginal[clock] = quadgk(t -> exp(step_likelihood(cs, clock, t)), cs.time, Inf)[1]
        end
    end
    total = sum(values(marginal))
    @assert abs(1.0 - total) < 1e-6
end
function execute(cs::ClockState, action::Fire)
    @assert action.clock ∈ keys(cs.clocks)
    @assert cs.clocks[action.clock].enabled
    if cs.check_sum
        core_matrix_sum(cs)
    end
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
    push!(cs.step_likelihood, steploglikelihood(cs.sampler, action.time, action.clock))
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
    if t1 < cs.clocks[fired_idx].te
        return -Inf
    end
    ll = 0.0
    @assert t0 < t1  # Assert invariant that time steps aren't simultaneous.
    for spec in values(cs.clocks)
        @assert spec.t0 <= t0  # Clock can't be enabled after last time point.
        !spec.enabled && continue
        te = spec.te
        if spec.t0 == te  # enabled at the zero-point of the distribution.
            if spec.clock == fired_idx
                @assert te < t1  # Can't fire unless pdf is nonzero.
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
                @assert te < t1  # Can't fire unless pdf is nonzero.
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
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_single_delay_clock = "Single delay clock" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.2),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_single_clock_gamma = "Single clock Gamma" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_single_delay_gamma = "Single delay clock Gamma" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), 0.2),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_single_soon_gamma = "Single soon clock Gamma" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), -0.2),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_multi_exp = "Multiple exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Enable(:c2, Exponential(1.2), 0.0),
        Enable(:c3, Exponential(0.3), 0.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_multi_nonexp = "Multiple non-exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Weibull(1.0), 0.0),
        Enable(:c2, Gamma(1.2), 0.0),
        Enable(:c3, Exponential(0.3), 0.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_multi__non_shift_exp = "Multiple exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Enable(:c2, Exponential(1.2), 0.2),
        Enable(:c3, Exponential(0.3), 1.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_multi_shift_nonexp = "Multiple non-exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Weibull(1.0), 0.0),
        Enable(:c2, Gamma(1.2), -0.2),
        Enable(:c3, Exponential(0.3), 0.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_multi__main_shift_exp = "Multiple exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.05),
        Enable(:c2, Exponential(1.2), 0.2),
        Enable(:c3, Exponential(0.3), 1.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_multi_main_shift_nonexp = "Multiple non-exponential" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Weibull(1.0), 0.05),
        Enable(:c2, Gamma(1.2), -0.2),
        Enable(:c3, Exponential(0.3), 0.0),
        Enable(:c4, Exponential(0.2), 0.0),
        Fire(:c1, 0.5),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_sequence = "Likely Sequence" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), 0.0),
        Fire(:c1, 0.5),
        Enable(:c2, Weibull(0.9), 0.0),
        Fire(:c2, 0.7),
        Enable(:c3, Exponential(1.2), 0.0),
        Fire(:c3, 1.2),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_resequence = "Likely Repeat Sequence" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), 0.0),
        Fire(:c1, 0.5),
        Enable(:c1, Weibull(0.9), 0.0),
        Fire(:c1, 0.7),
        Enable(:c1, Exponential(1.2), 0.0),
        Fire(:c1, 1.2),
    ]
    execute(cs, ss, actions)
end


@safetestset likely_sequence_reenable = "Likely Repeat Sequence Reenable" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Gamma(1.0), 0.0),
        Enable(:c2, Weibull(1.2), 0.2),
        Fire(:c1, 0.5),
        Enable(:c1, Weibull(0.9), 0.0),
        Enable(:c2, Gamma(0.9), 0.0),
        Fire(:c1, 0.7),
        Enable(:c1, Exponential(1.2), 0.0),
        Fire(:c2, 1.2),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_reenable_same_clock = "Re-enable same clock" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Fire(:c1, 0.5),
        # Re-enable the SAME clock with a DIFFERENT distribution
        Enable(:c1, Gamma(2.0, 0.5), 0.0),
        Fire(:c1, 1.2),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_disable_reenable = "Disable then re-enable" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        Enable(:c1, Exponential(1.0), 0.0),
        Enable(:c2, Exponential(2.0), 0.0),
        Fire(:c1, 0.3),
        # Disable c2 without it firing
        Disable(:c2),
        # Re-enable c2 with a different distribution
        Enable(:c2, Gamma(3.0, 0.2), 0.0),
        Fire(:c2, 0.8),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_time_dependent_offset = "Time-dependent with offset" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    actions = [
        # Enable at time 0 with offset, so distribution starts at 0.5
        Enable(:c1, Exponential(1.0), 0.5),
        # Enable another clock that can fire first
        Enable(:c2, Exponential(0.5), 0.0),
        Fire(:c2, 0.3),
        # Now re-enable c1 with a different offset relative to current time (0.3)
        Enable(:c1, Exponential(1.5), 0.2),  # Distribution starts at 0.3 + 0.2 = 0.5
        Fire(:c1, 0.7),
    ]
    execute(cs, ss, actions)
end

@safetestset likely_multiple_reenable_cycles = "Multiple re-enable cycles" begin
    using Distributions
    using ..LikelihoodHelper
    using CompetingClocks
    using Random
    rng = Xoshiro(2432343)
    K, T = (Symbol, Float64)
    builder = SamplerBuilder(K, T; sampler_spec=:firsttofire, trajectory_likelihood=true)
    sampler = SamplingContext(builder, rng)
    cs = ClockState()
    cs.check_sum = true
    ss = SampleState(sampler)
    # Simulate the pattern: enable transcribe, fire, re-enable transcribe
    actions = [
        Enable(:transcribe, Exponential(10.0), 0.0),
        Enable(:degrade, Gamma(4.0, 0.25), 0.0),
        Fire(:transcribe, 0.1),
        # After transcribe fires, re-enable it
        Enable(:transcribe, Exponential(10.0), 0.0),
        Fire(:degrade, 0.2),
        Fire(:transcribe, 0.25),
        Enable(:transcribe, Exponential(10.0), 0.0),
        Fire(:transcribe, 0.4),
    ]
    execute(cs, ss, actions)
end
