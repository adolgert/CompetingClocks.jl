# Gen.jl integration with CompetingClocks - Observation Likelihood Pattern
#
# This example demonstrates computing observation likelihoods for event data
# using CompetingClocks.pathloglikelihood combined with Gen's @factor macro.
#
# Statistical context: Given observed event times from a continuous-time system,
# we want to infer the parameters θ. The log-likelihood has the form:
#   log p(D|θ) = Σ log λ(t_i) - ∫ Λ(t) dt
# where λ is the hazard and Λ is the total hazard over all enabled clocks.
#
# CompetingClocks handles the hazard integration internally via pathloglikelihood().
#
# Run with: julia --project=examples examples/gen_observation.jl

println("Loading packages...")
using Random
using Distributions
using CompetingClocks
println("Packages loaded.")

# ===========================================================================
# 1. Data representation
# ===========================================================================

# Each observed event: `key` is :fail or :repair, `time` is event time.
const RelEvent = NamedTuple{(:key, :time), Tuple{Symbol, Float64}}

# Example toy dataset (normally you would read this from your data)
observed_events = RelEvent[
    (key = :fail,   time = 1.2),
    (key = :repair, time = 2.0),
    (key = :fail,   time = 4.5),
    (key = :repair, time = 5.1),
]
obs_end_time = 8.0  # last observation time T


# ===========================================================================
# 2. Reliability model state
# ===========================================================================

struct ReliabilityParams
    k_fail::Float64     # Weibull shape
    η_fail::Float64     # Weibull scale
    μ_repair::Float64   # Exponential repair rate (1/mean)
end

mutable struct ReliabilityModel
    up::Bool            # true if machine is working
    θ::ReliabilityParams
end

function init_model(θ::ReliabilityParams)
    ReliabilityModel(true, θ)  # start in the working state at t = 0
end


# ===========================================================================
# 3. Enabling clocks in CompetingClocks
# ===========================================================================

# Given the current state, enable whichever clock should be active
# Note: `now` parameter reserved for future use with time-dependent distributions
function enable_current_clock!(model::ReliabilityModel, sampler, now::Float64=0.0)

    if model.up
        # Time to next failure while machine is up
        dist = Weibull(model.θ.k_fail, model.θ.η_fail)
        enable!(sampler, :fail, dist)
    else
        # Time to next repair while machine is down
        # Exponential(μ) has mean μ, so rate = 1/μ
        dist = Exponential(model.θ.μ_repair)
        enable!(sampler, :repair, dist)
    end
end

# Advance the model by a known event at known time, updating state and clocks
function step_reliability!(model::ReliabilityModel,
                           sampler,
                           key::Symbol,
                           when::Float64)

    # Tell CompetingClocks that this event fired at "when"
    fire!(sampler, key, when)

    # Update physical state
    if key == :fail
        model.up = false
    elseif key == :repair
        model.up = true
    else
        error("Unknown event key: $key")
    end

    # Enable the next event appropriate for the new state
    enable_current_clock!(model, sampler, when)
end


# ===========================================================================
# 4. Pure Julia log-likelihood function
# ===========================================================================

"""
    reliability_path_loglikelihood(θ, events, end_time; rng=Xoshiro(1))

Compute log p(events | θ) for the reliability model over [0, end_time]
using CompetingClocks.pathloglikelihood.
"""
function reliability_path_loglikelihood(θ::ReliabilityParams,
                                        events::Vector{RelEvent},
                                        end_time::Float64;
                                        rng::AbstractRNG = Xoshiro(1))

    # Build a sampler with path-likelihood enabled.
    # `path_likelihood=true` means we can call `pathloglikelihood` later.
    sampler = SamplingContext(Symbol, Float64, rng;
                              path_likelihood = true)

    # Initialize model at time 0 and enable the first relevant clock
    model = init_model(θ)
    enable_current_clock!(model, sampler, 0.0)

    # Replay observed events in time order
    for evt in events
        key = evt.key
        t   = evt.time
        step_reliability!(model, sampler, key, t)
    end

    # Log-likelihood of whole path, including no further events before end_time
    log_prob = pathloglikelihood(sampler, end_time)

    return log_prob
end


# ===========================================================================
# Test Part 1: Basic likelihood computation
# ===========================================================================

println("=" ^ 60)
println("Part 1: Basic likelihood computation")
println("=" ^ 60)

# True parameters for testing
θ_true = ReliabilityParams(1.5, 2.0, 1.0)  # k=1.5, η=2.0, μ=1.0

println("Observed events:")
for (i, evt) in enumerate(observed_events)
    println("  $i: $(evt.key) at t=$(evt.time)")
end
println("Observation window: [0, $obs_end_time]")

# Compute log-likelihood at true parameters
ll_true = reliability_path_loglikelihood(θ_true, observed_events, obs_end_time)
println("\nLog-likelihood at θ_true = $ll_true")

# Compare with different parameters
θ_alt = ReliabilityParams(1.0, 3.0, 0.5)  # Different parameters
ll_alt = reliability_path_loglikelihood(θ_alt, observed_events, obs_end_time)
println("Log-likelihood at θ_alt  = $ll_alt")
println()


# ===========================================================================
# Part 2: Gen model using custom distribution for likelihood
# ===========================================================================

println("Loading Gen...")
using Gen
using Statistics: mean, std
println("Gen loaded.")

println("=" ^ 60)
println("Part 2: Gen model with likelihood factor")
println("=" ^ 60)

# Custom distribution that contributes a fixed log-likelihood to the trace.
# This acts like @factor() by adding an external likelihood term.
# The "sampled value" is just a dummy (nothing).
struct LikelihoodFactor <: Gen.Distribution{Nothing} end
const likelihood_factor = LikelihoodFactor()

function Gen.random(::LikelihoodFactor, logpdf_val::Float64)
    return nothing
end

function Gen.logpdf(::LikelihoodFactor, value::Nothing, logpdf_val::Float64)
    return logpdf_val
end

Gen.is_discrete(::LikelihoodFactor) = true
Gen.has_output_grad(::LikelihoodFactor) = false
Gen.has_argument_grads(::LikelihoodFactor) = (false,)

# A Gen generative function: sample θ, then add the CompetingClocks log-likelihood
@gen function reliability_gen_model(events::Vector{RelEvent},
                                    end_time::Float64)

    # Priors on log-parameters (unconstrained). Adjust as needed.
    log_k ~ normal(0.0, 0.5)   # shape > 0
    log_η ~ normal(0.0, 1.0)   # scale > 0
    log_μ ~ normal(0.0, 1.0)   # rate > 0

    k   = exp(log_k)
    η   = exp(log_η)
    μ   = exp(log_μ)

    θ = ReliabilityParams(k, η, μ)

    # External likelihood term from CompetingClocks
    loglike = reliability_path_loglikelihood(θ, events, end_time)

    # Add this as a factor to the Gen trace using our custom distribution
    {:likelihood} ~ likelihood_factor(loglike)

    # Return parameters so we can look at the posterior
    return θ
end


# ===========================================================================
# Test Part 2: Generate initial trace
# ===========================================================================

println("\nGenerating initial trace...")

# Arguments are (events, end_time)
(trace, weight) = generate(reliability_gen_model, (observed_events, obs_end_time))

println("Initial sampled params: ", get_retval(trace))
println("Initial log-joint (model + likelihood): ", get_score(trace))
println("Importance weight: $weight")
println()


# ===========================================================================
# Part 3: Simple Metropolis-Hastings inference
# ===========================================================================

println("=" ^ 60)
println("Part 3: Metropolis-Hastings inference")
println("=" ^ 60)

# Simple Gaussian drift MH for the three log-parameters
addrs = [:log_k, :log_η, :log_μ]

# Run MH (wrapped in function to avoid Julia scoping issues)
function run_mh(initial_trace, n_iters; drift_std=0.1)
    samples = Vector{ReliabilityParams}(undef, n_iters)
    scores = Vector{Float64}(undef, n_iters)
    current_trace = initial_trace
    accepted = 0

    for i in 1:n_iters
        # Cycle through parameters, proposing one at a time
        addr = addrs[(i - 1) % 3 + 1]

        # Get current value and propose new one
        current_val = current_trace[addr]
        proposed_val = current_val + randn() * drift_std

        # Create constraint with proposed value
        constraints = choicemap((addr, proposed_val))

        # Use update with constraints for deterministic proposal
        new_trace, weight, _, _ = update(current_trace, get_args(current_trace),
                                         (NoChange(), NoChange()), constraints)

        # Accept/reject (symmetric proposal, so weight is log(p_new/p_old))
        if log(rand()) < weight
            current_trace = new_trace
            accepted += 1
        end

        samples[i] = get_retval(current_trace)
        scores[i] = get_score(current_trace)
    end

    return samples, scores, accepted, current_trace
end

println("Running Metropolis-Hastings for 1000 iterations...")
n_iters = 1000
samples, scores, accepted, final_trace = run_mh(trace, n_iters)

println("Acceptance rate: $(accepted / n_iters)")
println()

# Compute posterior statistics (discard first half as burn-in)
burn_in = n_iters ÷ 2
post_samples = samples[burn_in+1:end]

k_samples = [s.k_fail for s in post_samples]
η_samples = [s.η_fail for s in post_samples]
μ_samples = [s.μ_repair for s in post_samples]

println("Posterior estimates (after burn-in):")
println("  k (Weibull shape): mean=$(round(mean(k_samples), digits=3)), std=$(round(std(k_samples), digits=3))")
println("  η (Weibull scale): mean=$(round(mean(η_samples), digits=3)), std=$(round(std(η_samples), digits=3))")
println("  μ (Repair mean):   mean=$(round(mean(μ_samples), digits=3)), std=$(round(std(μ_samples), digits=3))")
println()

println("True parameters for reference:")
println("  k = $(θ_true.k_fail), η = $(θ_true.η_fail), μ = $(θ_true.μ_repair)")
println()


# ===========================================================================
# Part 4: Generate synthetic data and verify inference
# ===========================================================================

println("=" ^ 60)
println("Part 4: Synthetic data generation and inference")
println("=" ^ 60)

"""
Generate synthetic events from the reliability model.
"""
function generate_synthetic_data(θ::ReliabilityParams, end_time::Float64;
                                 rng::AbstractRNG = Xoshiro(42))
    sampler = SamplingContext(Symbol, Float64, rng)
    model = init_model(θ)
    events = RelEvent[]

    enable_current_clock!(model, sampler, 0.0)

    when, which = next(sampler)
    while !isnothing(which) && when <= end_time
        fire!(sampler, which, when)
        push!(events, (key = which, time = when))

        # Update state
        if which == :fail
            model.up = false
        else
            model.up = true
        end
        enable_current_clock!(model, sampler, when)

        when, which = next(sampler)
    end

    return events
end

# Generate synthetic data with known parameters
θ_synth = ReliabilityParams(2.0, 3.0, 1.5)  # k=2.0, η=3.0, μ=1.5
synth_end_time = 50.0
synth_events = generate_synthetic_data(θ_synth, synth_end_time)

println("Generated $(length(synth_events)) synthetic events with:")
println("  k = $(θ_synth.k_fail), η = $(θ_synth.η_fail), μ = $(θ_synth.μ_repair)")
println()

# Run inference on synthetic data
println("Running MH on synthetic data (500 iterations)...")
(synth_trace, _) = generate(reliability_gen_model, (synth_events, synth_end_time))

n_synth_iters = 500
synth_samples, _, _, _ = run_mh(synth_trace, n_synth_iters)

# Posterior statistics
synth_post = synth_samples[n_synth_iters ÷ 2 + 1:end]
k_synth = [s.k_fail for s in synth_post]
η_synth = [s.η_fail for s in synth_post]
μ_synth = [s.μ_repair for s in synth_post]

println("\nPosterior estimates:")
println("  k: mean=$(round(mean(k_synth), digits=3)) (true: $(θ_synth.k_fail))")
println("  η: mean=$(round(mean(η_synth), digits=3)) (true: $(θ_synth.η_fail))")
println("  μ: mean=$(round(mean(μ_synth), digits=3)) (true: $(θ_synth.μ_repair))")
println()

println("=" ^ 60)
println("Example complete!")
println("=" ^ 60)
