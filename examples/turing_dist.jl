# Demonstrates using CompetingClocks.jl with Turing.jl for Bayesian inference
# on continuous-time stochastic process parameters.
#
# This example implements a birth-death model where:
# - Births occur at rate β * N(t) (Exponential distribution)
# - Deaths occur with Gamma-distributed lifetimes
#
# We observe a complete event path and infer the parameters (β, shape, scale)
# using MCMC via Turing.jl.

using Random
using Distributions
using Turing
using CompetingClocks

# ----------------------------
# 1. Birth–death model (CompetingClocks)
# ----------------------------

# Clock keys: (:birth, 0) for global birth, (:death, id) for individual deaths.
const ClockKey = Tuple{Symbol,Int}

struct BDParams
    β::Float64       # per-individual birth rate
    shape::Float64   # Gamma shape for lifetime
    scale::Float64   # Gamma scale for lifetime
end

mutable struct BDState
    population::Set{Int}  # current living individuals
    next_id::Int          # next ID to assign on birth
end

function BDState(N0::Int)
    pop = Set(1:N0)
    BDState(pop, N0 + 1)
end

"""
    initialize!(state, sampler, params)

Enable initial birth and death clocks at time 0.
"""
function initialize!(state::BDState, sampler, params::BDParams)
    N = length(state.population)
    N == 0 && return

    # Global birth clock: rate β * N ⇒ Exponential(scale = 1 / (β N))
    enable!(sampler, (:birth, 0),
            Exponential(1 / (params.β * N)))

    # One Gamma-distributed death clock per individual
    for id in state.population
        enable!(sampler, (:death, id),
                Gamma(params.shape, params.scale))
    end
end

"""
    handle_event!(state, sampler, evt, t, params)

Apply a birth or death event to state and update enabled clocks.
"""
function handle_event!(state::BDState,
                       sampler,
                       evt::ClockKey,
                       t::Float64,
                       params::BDParams)

    kind, idx = evt

    if kind === :birth
        new_id = state.next_id
        state.next_id += 1
        push!(state.population, new_id)

        # New individual's death clock
        enable!(sampler, (:death, new_id),
                Gamma(params.shape, params.scale))

    elseif kind === :death
        # fire!(...) already disabled (:death, idx)
        delete!(state.population, idx)

    else
        error("Unknown event kind $kind")
    end

    # Update global birth clock (depends on current N(t))
    if !isempty(state.population)
        enable!(sampler, (:birth, 0),
                Exponential(1 / (params.β * length(state.population))))
    else
        # No individuals left ⇒ no more births
        disable!(sampler, (:birth, 0))
    end
end

# ----------------------------
# 2. Path type and simulator
# ----------------------------

struct BDPath
    events::Vector{ClockKey}   # (:birth,0) or (:death,id)
    times::Vector{Float64}     # strictly increasing
    T_end::Float64             # observation horizon
end

BDPath(events::Vector{ClockKey}, times::Vector{Float64}) =
    BDPath(events, times, isempty(times) ? 0.0 : maximum(times))

Base.length(p::BDPath) = length(p.events)

# Required for Turing's model checking
Base.isnan(p::BDPath) = any(isnan, p.times) || isnan(p.T_end)

"""
    simulate_path(params; N0=10, T_end=10.0, rng=Xoshiro(42))

Simulate a birth–death path with CompetingClocks.
"""
function simulate_path(params::BDParams;
                       N0::Int      = 10,
                       T_end::Float64 = 10.0,
                       rng          = Random.Xoshiro(42))

    state   = BDState(N0)
    sampler = SamplingContext(ClockKey, Float64, rng)

    initialize!(state, sampler, params)

    events = ClockKey[]
    times  = Float64[]

    while true
        when, which = next(sampler)
        if isnothing(which) || when > T_end
            break
        end

        push!(events, which)
        push!(times,  when)

        fire!(sampler, which, when)
        handle_event!(state, sampler, which, when, params)
    end

    return BDPath(events, times, T_end)
end

# ----------------------------
# 3. Path likelihood as a distribution
# ----------------------------

"""
Distribution over whole event paths for the birth–death model.
The "sample" type is BDPath; parameters are β, shape, scale, and N0.
"""
struct BirthDeathPathDist <: Distributions.ContinuousMultivariateDistribution
    β::Float64
    shape::Float64
    scale::Float64
    N0::Int
end

function Distributions.insupport(d::BirthDeathPathDist, path::BDPath)
    # Very minimal check; you can add monotonicity checks, etc.
    length(path.events) == length(path.times)
end

"""
    Distributions.logpdf(d::BirthDeathPathDist, path::BDPath)

Replay `path` into a SamplingContext with `path_likelihood=true`
and return the CompetingClocks path log-likelihood.
"""
function Distributions.logpdf(d::BirthDeathPathDist, path::BDPath)
    @assert Distributions.insupport(d, path)

    params = BDParams(d.β, d.shape, d.scale)

    # Fresh sampler with path likelihood enabled.
    # RNG is irrelevant since we never call `next`, only `fire!`.
    rng     = Random.Xoshiro(0)
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)

    state = BDState(d.N0)
    initialize!(state, sampler, params)

    # Replay the observed path deterministically
    for (evt, t) in zip(path.events, path.times)
        fire!(sampler, evt, t)
        handle_event!(state, sampler, evt, t, params)
    end

    # Full path log-likelihood, including survival to T_end
    # pathloglikelihood returns a scalar when using single distribution
    ll = pathloglikelihood(sampler, path.T_end)
    return ll isa AbstractVector ? ll[1] : ll
end

# Required for ContinuousMultivariateDistribution
Distributions.length(d::BirthDeathPathDist) = 1  # Not really used, but required

# Required for Turing - loglikelihood is the same as logpdf for a single observation
function Distributions.loglikelihood(d::BirthDeathPathDist, path::BDPath)
    Distributions.logpdf(d, path)
end

# ----------------------------
# 4. Turing model
# ----------------------------

@model function birthdeath_model(path::BDPath, N0::Int)
    # Weakly informative priors around reasonable values
    # True values: β=0.3, shape=2.0, scale=0.5
    β     ~ LogNormal(log(0.3), 0.5)
    shape ~ LogNormal(log(2.0), 0.5)
    scale ~ LogNormal(log(0.5), 0.5)

    # Likelihood: the whole path is drawn from the CT birth–death model
    path  ~ BirthDeathPathDist(β, shape, scale, N0)
end

# ----------------------------
# 5. Synthetic data and inference
# ----------------------------

function main()
    println("Birth-Death Process Inference with CompetingClocks + Turing")
    println("=" ^ 60)

    # "True" parameters for data generation
    # Note: For stable simulation, need birth rate < death rate.
    # With Gamma(shape, scale), mean lifetime = shape * scale.
    # Death rate per individual ≈ 1/(mean lifetime).
    # For stability: β < 1/(shape * scale)
    # Here: β=0.3, Gamma(2, 0.5), mean=1.0, death_rate≈1.0
    # So birth rate (0.3) < death rate (1.0) → population declines
    true_params = BDParams(0.3, 2.0, 0.5)
    N0          = 20
    T_end       = 10.0
    rng_data    = Random.Xoshiro(1234)

    println("\nTrue parameters:")
    println("  β     = $(true_params.β)")
    println("  shape = $(true_params.shape)")
    println("  scale = $(true_params.scale)")

    println("\nSimulating observed path...")
    observed_path = simulate_path(true_params; N0=N0, T_end=T_end, rng=rng_data)
    println("  Observed $(length(observed_path)) events over [0, $T_end]")

    # Count event types
    births = count(e -> e[1] === :birth, observed_path.events)
    deaths = count(e -> e[1] === :death, observed_path.events)
    println("  Births: $births, Deaths: $deaths")

    # Test the likelihood calculation
    println("\nTesting likelihood calculation...")
    test_ll = Distributions.logpdf(
        BirthDeathPathDist(true_params.β, true_params.shape, true_params.scale, N0),
        observed_path
    )
    println("  Log-likelihood at true parameters: $(round(test_ll, digits=2))")

    # Build Turing model
    println("\nRunning MCMC inference...")
    model = birthdeath_model(observed_path, N0)

    # Run MH sampler (simpler than NUTS for this example)
    rng_mcmc = Random.Xoshiro(2024)
    chain = sample(model, MH(), 500; rng=rng_mcmc, progress=true)

    println("\n" * "=" ^ 60)
    println("MCMC Results:")
    println(chain)

    # Summary statistics
    println("\nParameter recovery:")
    for (param, true_val) in [("β", true_params.β), ("shape", true_params.shape), ("scale", true_params.scale)]
        samples = chain[Symbol(param)].data
        est_mean = mean(samples)
        est_std = std(samples)
        println("  $param: true=$(round(true_val, digits=3)), " *
                "estimated=$(round(est_mean, digits=3)) ± $(round(est_std, digits=3))")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
