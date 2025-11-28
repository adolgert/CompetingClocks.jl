# Gen.jl integration with CompetingClocks using a custom Generative Function
#
# This example demonstrates how to wrap a CompetingClocks simulator as a
# Gen custom generative function, enabling it to be used as a building block
# in larger probabilistic models.
#
# Run with: julia --project=examples examples/gen_generative.jl

println("Loading packages...")
using Random
using Distributions
using CompetingClocks
println("Packages loaded.")

# ===========================================================================
# Part 1: Birth-death process simulation with path likelihood
# ===========================================================================

# Clock key layout:
#   (:birth, 0) for birth events
#   (:death, i) for the death of individual i
const ClockKey = Tuple{Symbol,Int}

# One event in a path
struct BDEvent
    time::Float64
    key::ClockKey
end

const EventPath = Vector{BDEvent}

"""
    simulate_bd_path(birth_rate, K, death_shape, death_scale, init_pop, tmax, rng)

Run one CompetingClocks simulation of a density-dependent birth-death process.
Returns (path, final_population, log_path_prob).

The birth rate follows a logistic form:
    birth_rate_total = birth_rate * N * max(0, 1 - N/K)

This ensures stable population dynamics around carrying capacity K.
"""
function simulate_bd_path(
    birth_rate::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    init_pop::Int,
    tmax::Float64,
    rng::AbstractRNG,
)
    # Initial state: individuals are numbered 1:init_pop
    population = Set(1:init_pop)
    next_id = init_pop + 1

    # CompetingClocks sampler with path likelihood enabled
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)

    # Density-dependent birth rate: logistic model
    function birth_rate_total(N::Int)
        N <= 0 && return 0.0
        return birth_rate * N * max(0.0, 1.0 - N / K)
    end

    # Enable initial birth and death clocks
    rate = birth_rate_total(length(population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Simulate until horizon tmax or extinction
    path = BDEvent[]
    when, which = next(sampler)

    while !isnothing(which) && when <= tmax && !isempty(population)
        # Advance sampler
        fire!(sampler, which, when)

        # Record event
        push!(path, BDEvent(when, which))

        if which[1] == :birth
            # New individual
            new_id = next_id
            next_id += 1
            push!(population, new_id)
            # Enable its death clock
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif which[1] == :death
            delete!(population, which[2])
        else
            error("Unexpected clock key $(which)")
        end

        # Update the birth clock with density-dependent rate
        rate = birth_rate_total(length(population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end
        # If rate == 0 (at or above K), birth clock stays disabled

        when, which = next(sampler)
    end

    # Log-likelihood of the whole path, including survival up to tmax
    logp = pathloglikelihood(sampler, tmax)

    return path, length(population), logp
end

# ===========================================================================
# Test Part 1: Basic simulation
# ===========================================================================

println("=" ^ 60)
println("Part 1: Birth-death simulation with path likelihood")
println("=" ^ 60)

rng = Xoshiro(12345)
birth_rate = 2.0
K = 50.0  # Carrying capacity
death_shape = 2.0
death_scale = 2.0  # mean death time = shape * scale = 4.0
init_pop = 10
tmax = 10.0

path, final_pop, logp = simulate_bd_path(
    birth_rate, K, death_shape, death_scale, init_pop, tmax, rng
)

println("Simulated $(length(path)) events in [0, $tmax]")
println("Final population: $final_pop")
println("Path log-likelihood: $logp")
println()

# ===========================================================================
# Part 2: Custom Gen generative function
# ===========================================================================

println("Loading Gen...")
using Gen
println("Gen loaded.")

println("=" ^ 60)
println("Part 2: Custom Gen generative function")
println("=" ^ 60)

# Trace type: stores the simulation results
struct BDTrace <: Gen.Trace
    args::Tuple                # (birth_rate, K, death_shape, death_scale, init_pop, tmax)
    choices::Gen.ChoiceMap     # :path => Vector{BDEvent}
    retval::Int                # final population size
    logp::Float64              # log p(path | args)
    gen_fn::Any                # back-pointer to the generative function
end

# Generative function type: returns Int (final population), uses BDTrace as trace
struct BDPathGF{R<:AbstractRNG} <: Gen.GenerativeFunction{Int,BDTrace}
    rng::R
end

BDPathGF(rng::AbstractRNG) = BDPathGF{typeof(rng)}(rng)
BDPathGF() = BDPathGF(Random.default_rng())

# Implement the GFI methods required for a custom GF

function Gen.simulate(gen_fn::BDPathGF, args::Tuple)
    @assert length(args) == 6 "Expected 6 arguments: (birth_rate, K, death_shape, death_scale, init_pop, tmax)"

    birth_rate, K, death_shape, death_scale, init_pop, tmax = args

    path, final_pop, logp = simulate_bd_path(
        birth_rate, K, death_shape, death_scale, init_pop, tmax, gen_fn.rng
    )

    choices = Gen.choicemap((:path, path))
    return BDTrace(args, choices, final_pop, logp, gen_fn)
end

# Default proposal q = p, so generate just calls simulate and returns weight 0
function Gen.generate(gen_fn::BDPathGF, args::Tuple, constraints::Gen.ChoiceMap)
    trace = Gen.simulate(gen_fn, args)
    return trace, 0.0
end

Gen.generate(gen_fn::BDPathGF, args::Tuple) =
    Gen.generate(gen_fn, args, Gen.choicemap())

# Accessors for the trace
Gen.get_args(trace::BDTrace) = trace.args
Gen.get_retval(trace::BDTrace) = trace.retval
Gen.get_choices(trace::BDTrace) = trace.choices
Gen.get_score(trace::BDTrace) = trace.logp
Gen.get_gen_fn(trace::BDTrace) = trace.gen_fn

Base.getindex(trace::BDTrace, addr) = getindex(trace.choices, addr)

# No gradients in this example
Gen.has_argument_grads(::BDPathGF) = (false, false, false, false, false, false)
Gen.accepts_output_grad(::BDPathGF) = false

# Simple project implementation: return full log-probability
Gen.project(trace::BDTrace, ::Gen.Selection) = trace.logp

# ===========================================================================
# Test Part 2: Standalone generative function
# ===========================================================================

bd_gf = BDPathGF(Xoshiro(42))

trace = Gen.simulate(bd_gf, (birth_rate, K, death_shape, death_scale, init_pop, tmax))
final_pop_gf = Gen.get_retval(trace)
path_gf = trace[:path]
logp_gf = Gen.get_score(trace)

println("Standalone generative function test:")
println("  Final population: $final_pop_gf")
println("  Number of events: $(length(path_gf))")
println("  Log p(path | params): $logp_gf")
println()

# ===========================================================================
# Part 3: Using BDPathGF inside a Gen model
# ===========================================================================

println("=" ^ 60)
println("Part 3: Using BDPathGF inside a Gen model")
println("=" ^ 60)

# Create a global instance of the generative function
const global_bd_gf = BDPathGF(Xoshiro(123))

"""
A Gen model that:
- Takes carrying capacity K and horizon tmax as fixed parameters
- Puts a prior on the birth rate
- Runs the CompetingClocks simulation as a sub-generative-function
- Observes noisy measurement of final population
"""
@gen function bd_inference_model(K::Float64, tmax::Float64)
    # Prior on the birth rate (lognormal via exp of normal)
    log_birth_rate = @trace(normal(0.0, 0.5), :log_birth_rate)
    birth_rate = exp(log_birth_rate)

    # Fixed death parameters and initial population
    death_shape = 2.0
    death_scale = 2.0
    init_pop = 10

    # Draw a whole birth-death path via CompetingClocks.
    # Return value is the final population size at tmax.
    final_pop = @trace(
        global_bd_gf(birth_rate, K, death_shape, death_scale, init_pop, tmax),
        :simulation
    )

    # Noisy observation of the final population
    y = @trace(normal(float(final_pop), 2.0), :y)

    return (birth_rate=birth_rate, log_birth_rate=log_birth_rate, final_pop=final_pop, y=y)
end

# ===========================================================================
# Test Part 3: Forward simulation from the model
# ===========================================================================

println("Forward simulation from bd_inference_model:")
trace_model = Gen.simulate(bd_inference_model, (K, tmax))
retval = Gen.get_retval(trace_model)

println("  Sampled birth_rate: $(retval.birth_rate)")
println("  Final population: $(retval.final_pop)")
println("  Noisy observation y: $(retval.y)")
println("  Model log-probability: $(Gen.get_score(trace_model))")
println()

# ===========================================================================
# Part 4: Inference - conditioning on observed data
# ===========================================================================

println("=" ^ 60)
println("Part 4: Inference with importance sampling")
println("=" ^ 60)

# Generate synthetic observed data
println("Generating synthetic observation...")
true_birth_rate = 1.5
true_bd_gf = BDPathGF(Xoshiro(999))
_, true_final_pop, _ = simulate_bd_path(
    true_birth_rate, K, death_shape, death_scale, init_pop, tmax, Xoshiro(999)
)
# Add noise to get observation
observed_y = true_final_pop + randn() * 2.0

println("  True birth rate: $true_birth_rate")
println("  True final population: $true_final_pop")
println("  Observed y: $observed_y")
println()

# Condition on the observed y
observations = Gen.choicemap((:y, observed_y))

# Run importance sampling
println("Running importance sampling (100 samples)...")
n_samples = 100
traces, log_weights, lml_est = Gen.importance_sampling(
    bd_inference_model, (K, tmax), observations, n_samples
)

# Analyze posterior
posterior_birth_rates = [exp(tr[:log_birth_rate]) for tr in traces]
posterior_final_pops = [Gen.get_retval(tr).final_pop for tr in traces]

# Compute weighted statistics
weights = exp.(log_weights .- maximum(log_weights))
weights ./= sum(weights)

mean_birth_rate = sum(weights .* posterior_birth_rates)
mean_final_pop = sum(weights .* posterior_final_pops)

println("Posterior summary (importance-weighted):")
println("  Mean birth rate: $(round(mean_birth_rate, digits=3)) (true: $true_birth_rate)")
println("  Mean final pop: $(round(mean_final_pop, digits=1)) (observed: $(round(observed_y, digits=1)))")
println("  Log marginal likelihood estimate: $(round(lml_est, digits=2))")
println()

# Access choices from the trace
# Note: Nested custom GF choices require using Gen.get_choices and get_submap
example_trace = traces[1]
choices = Gen.get_choices(example_trace)
sim_choices = Gen.get_submap(choices, :simulation)
if !isempty(sim_choices)
    example_path = sim_choices[:path]
    println("Example path from first sample has $(length(example_path)) events")
else
    # With custom GFs, internal choices are stored in the sub-trace
    # The path log-likelihood is included in the model score
    println("(Path stored internally in custom generative function)")
end
println()

println("=" ^ 60)
println("Example complete!")
println("=" ^ 60)
