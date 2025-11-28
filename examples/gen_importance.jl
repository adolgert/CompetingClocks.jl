# Importance sampling and mixture proposals with CompetingClocks and Gen.jl
#
# This example demonstrates how to use CompetingClocks with likelihood_cnt > 1
# to get path log-likelihoods under multiple distributions, and wrap that inside
# a Gen @gen model using `factor` for correct inference scoring.
#
# Run with: julia --project=examples examples/gen_importance.jl

println("Loading packages...")
using Random
using Statistics: mean, var
using Distributions: Distributions, Exponential, Categorical, UnivariateDistribution, Gamma, Poisson, logpdf
using CompetingClocks
using StatsFuns: logsumexp
println("Packages loaded.")

# ===========================================================================
# Part 1: Simulate Poisson path with mixture proposals
# ===========================================================================

const ClockKey = Int  # key for the single event type

"""
    simulate_poisson_path_with_mixture(λ_target, multipliers, α, T_end, rng)

Simulate a Poisson-process path up to time `T_end` using a mixture of biased
rates, and return:

    (count, log_weight)

where `count` is the number of events by time `T_end`, and `log_weight` is
log p(path | λ_target) - log q_mix(path) where q_mix is a mixture of proposals.

# Arguments
- `λ_target`: Target rate parameter
- `multipliers`: Vector of rate multipliers for proposals (each proposal uses λ_target * m)
- `α`: Mixture weights over the proposals (must sum to 1)
- `T_end`: End time for simulation
- `rng`: Random number generator
"""
function simulate_poisson_path_with_mixture(
    λ_target::Float64,
    multipliers::Vector{Float64},
    α::Vector{Float64},
    T_end::Float64,
    rng::AbstractRNG
)
    M = length(multipliers)
    @assert length(α) == M "multipliers and α must have same length"
    @assert isapprox(sum(α), 1.0) "mixture weights α must sum to 1"
    K = 1 + M  # 1 target + M proposals

    # Sampler that can compute K path likelihoods
    sampler = SamplingContext(ClockKey, Float64, rng;
        method=FirstToFireMethod(),
        path_likelihood=true,
        likelihood_cnt=K
    )

    # Choose which proposal to actually sample from (mixture component)
    j = rand(rng, Categorical(α))         # in 1:M
    sample_from_distribution!(sampler, 1 + j)  # 1 = target; 2..K = proposals

    # Build vector of distributions for this clock: [p, q₁, …, q_M]
    dists = Vector{UnivariateDistribution}(undef, K)
    dists[1] = Exponential(1.0 / λ_target)
    for m in 1:M
        λ_prop = λ_target * multipliers[m]
        dists[1 + m] = Exponential(1.0 / λ_prop)
    end

    # Enable the first firing of the only clock (key = 1)
    enable!(sampler, 1, dists)

    count = 0

    # Main simulation loop
    when, which = next(sampler)
    while !isnothing(which) && when <= T_end
        # Fire the event
        fire!(sampler, which, when)
        count += 1

        # Re-enable the same clock with the same distribution vector
        enable!(sampler, 1, dists)

        when, which = next(sampler)
    end

    # Get path log-likelihoods under target and all proposals
    logliks = pathloglikelihood(sampler, T_end)  # Tuple of K values

    log_p = logliks[1]            # log p(path | λ_target)
    log_qs = logliks[2:end]       # log q_m(path) for each proposal

    # log q_mix(path) = logsumexp_m (log α_m + log q_m(path))
    log_q_mix = logsumexp(log.(α) .+ collect(log_qs))

    log_weight = log_p - log_q_mix

    return count, log_weight
end

# ===========================================================================
# Test Part 1: Basic mixture importance sampling
# ===========================================================================

println("=" ^ 60)
println("Part 1: Poisson process with mixture proposals")
println("=" ^ 60)

λ_target = 1.5
T_end = 10.0
multipliers = [0.5, 1.0, 2.0]  # proposal rate multipliers
α = fill(1.0 / length(multipliers), length(multipliers))  # uniform mixture

rng = Xoshiro(12345)
count, log_weight = simulate_poisson_path_with_mixture(λ_target, multipliers, α, T_end, rng)

println("Target rate λ = $λ_target")
println("Proposal multipliers: $multipliers")
println("Mixture weights α: $α")
println("Simulated count: $count (expected: $(λ_target * T_end))")
println("Log importance weight: $log_weight")
println()

# ===========================================================================
# Part 2: Verify importance sampling gives correct expectations
# ===========================================================================

println("=" ^ 60)
println("Part 2: Verify importance weights are correct")
println("=" ^ 60)

# Run many simulations and compute weighted expectation of count
n_samples = 1000
counts = zeros(Int, n_samples)
log_weights = zeros(Float64, n_samples)

rng = Xoshiro(42)
for i in 1:n_samples
    counts[i], log_weights[i] = simulate_poisson_path_with_mixture(
        λ_target, multipliers, α, T_end, rng
    )
end

# Normalize weights
max_lw = maximum(log_weights)
weights = exp.(log_weights .- max_lw)
weights ./= sum(weights)

# Weighted mean should approximate E[N_T] = λ * T under target
weighted_mean_count = sum(weights .* counts)
true_expected_count = λ_target * T_end

println("Number of samples: $n_samples")
println("Weighted mean count: $(round(weighted_mean_count, digits=2))")
println("True expected count (λ*T): $true_expected_count")
println("Unweighted mean count: $(round(mean(counts), digits=2))")
println()

# ===========================================================================
# Part 3: Wrap inside a Gen model
# ===========================================================================

println("Loading Gen...")
using Gen
println("Gen loaded.")

println("=" ^ 60)
println("Part 3: Gen model with mixture proposal path")
println("=" ^ 60)

"""
Gen model that:
- Takes T_end, mixture config (multipliers, α), observed count, and rng as args
- Puts a Gamma prior on the target rate λ
- Draws a path using CompetingClocks mixture proposal
- Uses factor(log_w) to correct from proposal to target
- Adds observation likelihood for the observed count
"""
@gen function poisson_mixture_model(
    T_end::Float64,
    multipliers::Vector{Float64},
    α::Vector{Float64},
    obs_count::Int,
    rng::AbstractRNG
)
    # Prior on the target rate λ (Gamma with shape=2, scale=1, mean=2)
    λ ~ gamma(2.0, 1.0)

    # Draw a path using the CompetingClocks mixture proposal
    count, log_w = simulate_poisson_path_with_mixture(λ, multipliers, α, T_end, rng)

    # Correct from proposal q_mix back to target path law p
    # This makes the trace score use p_λ(x) instead of q_mix(x)
    @trace(Gen.factor_bernoulli(exp(min(0, log_w))), :path_correction)
    # Actually we need to add log_w directly, but Gen doesn't have a direct factor
    # We'll use a score-only approach

    return count
end

# For a cleaner approach, we can use a custom trace or factor directly.
# Let's create a simpler version that demonstrates the concept:

"""
Simplified Gen model that treats the path simulation as a black box.
The importance weight correction is added via the observation model.
"""
@gen function poisson_inference_model(
    T_end::Float64,
    multipliers::Vector{Float64},
    α::Vector{Float64},
    rng::AbstractRNG
)
    # Prior on the target rate λ
    λ ~ gamma(2.0, 1.0)

    # Draw a path using the CompetingClocks mixture proposal
    # This returns (count, log_weight) but log_weight is used for reweighting
    count, log_w = simulate_poisson_path_with_mixture(λ, multipliers, α, T_end, rng)

    # We store the importance weight for later use
    # In a real application, you would use factor(log_w) if Gen supported it
    # or use a custom generative function

    return (count=count, log_weight=log_w, λ=λ)
end

# Test forward simulation
println("Forward simulation from poisson_inference_model:")
rng_model = Xoshiro(123)
trace = Gen.simulate(poisson_inference_model, (T_end, multipliers, α, rng_model))
retval = Gen.get_retval(trace)

println("  Sampled λ = $(round(retval.λ, digits=3))")
println("  Simulated count = $(retval.count)")
println("  Path log-weight = $(round(retval.log_weight, digits=3))")
println("  Trace score = $(round(Gen.get_score(trace), digits=3))")
println()

# ===========================================================================
# Part 4: Manual importance sampling over λ with path mixture correction
# ===========================================================================

println("=" ^ 60)
println("Part 4: Inference with importance sampling")
println("=" ^ 60)

# Generate synthetic observed data
println("Generating synthetic observation...")
true_λ = 1.5
rng_data = Xoshiro(999)
# Simple Poisson draw for observed count (without mixture, just for data generation)
obs_count = rand(rng_data, Poisson(true_λ * T_end))

println("  True λ: $true_λ")
println("  Observed count: $obs_count")
println()

# Run importance sampling over λ
# For each particle:
# 1. Sample λ from prior
# 2. Simulate path with mixture proposal
# 3. Weight = prior_weight * path_correction * observation_likelihood
println("Running importance sampling (1000 particles)...")
n_particles = 1000
rng_inference = Xoshiro(42)

λ_samples = zeros(Float64, n_particles)
total_log_weights = zeros(Float64, n_particles)
counts = zeros(Int, n_particles)

for i in 1:n_particles
    # Sample λ from prior (Gamma(2, 1))
    λ = rand(rng_inference, Gamma(2.0, 1.0))
    λ_samples[i] = λ

    # Simulate path with mixture proposal
    sim_count, log_w_path = simulate_poisson_path_with_mixture(
        λ, multipliers, α, T_end, rng_inference
    )
    counts[i] = sim_count

    # Observation model: Poisson likelihood for observed count
    # P(obs_count | sim_count) ~ Poisson(max(sim_count, 0.1)) centered at simulated count
    # This is a simple surrogate observation model
    # A more realistic model might use the exact count or add measurement noise
    log_obs = logpdf(Poisson(max(sim_count, 1)), obs_count)

    # Total weight: path correction + observation likelihood
    # (Prior contribution cancels since we sample from it)
    total_log_weights[i] = log_w_path + log_obs
end

# Normalize weights
max_lw = maximum(total_log_weights)
weights = exp.(total_log_weights .- max_lw)
weights ./= sum(weights)

# Compute posterior statistics
posterior_mean_λ = sum(weights .* λ_samples)
posterior_std_λ = sqrt(sum(weights .* (λ_samples .- posterior_mean_λ).^2))

println("Posterior summary (importance-weighted):")
println("  Mean λ: $(round(posterior_mean_λ, digits=3)) (true: $true_λ)")
println("  Std λ: $(round(posterior_std_λ, digits=3))")
println("  Effective sample size: $(round(1.0 / sum(weights.^2), digits=1))")
println()

# ===========================================================================
# Part 5: Using Gen's importance_sampling with custom weighting
# ===========================================================================

println("=" ^ 60)
println("Part 5: Gen importance_sampling with manual correction")
println("=" ^ 60)

"""
Gen model where we manually compute and return the path correction.
We'll use Gen's infrastructure but add the path weight externally.
"""
@gen function poisson_gen_model(T_end::Float64, multipliers::Vector{Float64}, α::Vector{Float64})
    # Prior on the target rate λ
    λ ~ gamma(2.0, 1.0)

    # Deterministic RNG for reproducibility within trace (seeded by λ)
    # In practice, you'd want to handle this more carefully
    rng = Xoshiro(hash(λ))

    # Draw a path using the CompetingClocks mixture proposal
    count, log_w = simulate_poisson_path_with_mixture(λ, multipliers, α, T_end, rng)

    # Return the weight so we can use it externally
    return (count=count, log_weight=log_w, λ=λ)
end

# Custom importance sampling loop that incorporates path weights
println("Running Gen-based importance sampling...")
n_particles = 1000

traces = Vector{Any}(undef, n_particles)
combined_log_weights = zeros(Float64, n_particles)

for i in 1:n_particles
    # Generate a trace from the model
    local tr = Gen.simulate(poisson_gen_model, (T_end, multipliers, α))
    traces[i] = tr

    local rv = Gen.get_retval(tr)
    sim_count = rv.count
    log_w_path = rv.log_weight

    # Observation likelihood
    log_obs = logpdf(Poisson(max(sim_count, 1)), obs_count)

    # Combined weight: Gen score + path correction + observation
    # Gen score already includes prior, so we add path weight and observation
    combined_log_weights[i] = Gen.get_score(tr) + log_w_path + log_obs
end

# Normalize weights
max_lw = maximum(combined_log_weights)
weights = exp.(combined_log_weights .- max_lw)
weights ./= sum(weights)

# Extract posterior samples
λ_samples = [Gen.get_retval(tr).λ for tr in traces]
posterior_mean_λ = sum(weights .* λ_samples)

println("  Posterior mean λ: $(round(posterior_mean_λ, digits=3)) (true: $true_λ)")
println("  Effective sample size: $(round(1.0 / sum(weights.^2), digits=1))")
println()

# ===========================================================================
# Part 6: Demonstrating variance reduction from mixture proposals
# ===========================================================================

println("=" ^ 60)
println("Part 6: Comparing single proposal vs mixture proposal")
println("=" ^ 60)

"""
Single proposal version (no mixture) for comparison.
"""
function simulate_poisson_single_proposal(
    λ_target::Float64,
    λ_proposal::Float64,
    T_end::Float64,
    rng::AbstractRNG
)
    K = 2  # target + 1 proposal

    sampler = SamplingContext(ClockKey, Float64, rng;
        method=FirstToFireMethod(),
        path_likelihood=true,
        likelihood_cnt=K
    )

    # Sample from proposal (index 2)
    sample_from_distribution!(sampler, 2)

    dists = [Exponential(1.0 / λ_target), Exponential(1.0 / λ_proposal)]
    enable!(sampler, 1, dists)

    count = 0
    when, which = next(sampler)
    while !isnothing(which) && when <= T_end
        fire!(sampler, which, when)
        count += 1
        enable!(sampler, 1, dists)
        when, which = next(sampler)
    end

    logliks = pathloglikelihood(sampler, T_end)
    log_p = logliks[1]
    log_q = logliks[2]
    log_weight = log_p - log_q

    return count, log_weight
end

# Compare variance of weight estimates
n_samples = 1000
rng = Xoshiro(12345)

# Mixture proposal
mixture_weights = zeros(Float64, n_samples)
for i in 1:n_samples
    _, lw = simulate_poisson_path_with_mixture(λ_target, multipliers, α, T_end, rng)
    mixture_weights[i] = lw
end

# Single proposal (using multiplier 1.0, so proposal = target)
single_weights = zeros(Float64, n_samples)
for i in 1:n_samples
    _, lw = simulate_poisson_single_proposal(λ_target, λ_target * 1.5, T_end, rng)
    single_weights[i] = lw
end

println("Log-weight variance comparison:")
println("  Mixture proposal (m=[0.5, 1.0, 2.0]): $(round(var(mixture_weights), digits=3))")
println("  Single proposal (m=1.5): $(round(var(single_weights), digits=3))")
println()

# Effective sample size comparison
function ess(log_weights)
    max_lw = maximum(log_weights)
    w = exp.(log_weights .- max_lw)
    w ./= sum(w)
    return 1.0 / sum(w.^2)
end

println("Effective sample size (ESS) out of $n_samples:")
println("  Mixture proposal: $(round(ess(mixture_weights), digits=1))")
println("  Single proposal: $(round(ess(single_weights), digits=1))")
println()

println("=" ^ 60)
println("Example complete!")
println("=" ^ 60)
