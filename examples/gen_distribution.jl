# Gen.jl integration with CompetingClocks
#
# This example demonstrates how to treat a CompetingClocks trajectory as a
# single Gen random choice, using pathloglikelihood for the log-density.
#
# Run with: julia --project=examples examples/gen_distribution.jl

#  Key differences from the documentation

#   1. SamplerBuilder vs SamplingContext: The doc uses SamplerBuilder + SamplingContext(builder, rng):
#   # Doc version:
#   builder = SamplerBuilder(ClockKey, Float64; method=FirstToFireMethod(), path_likelihood=true)
#   sampler = SamplingContext(builder, rng)
#   1. I used the simpler direct constructor:
#   # My version:
#   sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)
#   1. Both work - the direct constructor is a convenience wrapper.
#   2. Function organization: The doc splits helpers into separate functions (enable_initial_clocks!, handle_event!). I inlined them for compactness to
#   reduce JIT compilation time (which was causing OOM).
#   3. Parameters: I used smaller values (N0=3, t_max=5.0) vs doc's (N0=10, t_max=10.0) to reduce memory during compilation.
#   4. Constant name: Doc uses bd_path, I used bd_path_dist (trivial difference).

println("Loading packages...")
using Random
using Distributions
using CompetingClocks
println("Packages loaded.")

# ===========================================================================
# Part 1: Basic birth-death process simulation
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

# State: set of live individuals and a counter for the next ID
mutable struct BDState
    population::Set{Int}
    next_id::Int
end

BDState(N0::Int) = BDState(Set(1:N0), N0 + 1)

"""
Run a birth-death simulation with density-dependent (logistic) birth rate.

The birth rate follows a logistic form:
    birth_rate = λ_birth * N * max(0, 1 - N/K)

This ensures the population fluctuates around carrying capacity K rather than
exploding or dying out.

Parameters:
- rng: Random number generator
- t_max: Maximum simulation time
- λ_birth: Per-capita birth rate at low density
- K: Carrying capacity (population where birth rate → 0)
- death_shape, death_scale: Gamma distribution parameters for death times
- N0: Initial population size
"""
function simulate_bd(
    rng::AbstractRNG,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)
    # Create sampler with path-likelihood enabled
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)
    state = BDState(N0)
    path = EventPath()

    # Density-dependent birth rate: logistic model
    function birth_rate(N::Int)
        N <= 0 && return 0.0
        return λ_birth * N * max(0.0, 1.0 - N / K)
    end

    # Enable initial clocks
    rate = birth_rate(length(state.population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in state.population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Simulation loop
    when, which = next(sampler)
    while !isnothing(which) && when <= t_max
        fire!(sampler, which, when)
        push!(path, BDEvent(when, which))

        # Handle event
        if which[1] == :birth
            new_id = state.next_id
            state.next_id += 1
            push!(state.population, new_id)
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif which[1] == :death
            delete!(state.population, which[2])
        end

        # Refresh birth clock with density-dependent rate
        rate = birth_rate(length(state.population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end
        # If rate == 0 (at or above K), birth clock stays disabled

        when, which = next(sampler)
    end

    return path, sampler
end

# ===========================================================================
# Test Part 1: Basic simulation
# ===========================================================================

println("=" ^ 60)
println("Part 1: Basic birth-death simulation")
println("=" ^ 60)

# Use small population to reduce compilation time
rng = Xoshiro(12345)
t_max = 10.0
λ_birth = 2.0
K = 50.0  # Carrying capacity - population fluctuates around this
death_shape = 2.0
death_scale = 2.0  # mean death time = shape * scale = 4.0
N0 = 10

path, sampler = simulate_bd(rng, t_max, λ_birth, K, death_shape, death_scale, N0)

println("Simulated $(length(path)) events in [0, $t_max]")
println("First few events:")
for (i, ev) in enumerate(first(path, 5))
    println("  $i: time=$(round(ev.time, digits=4)), key=$(ev.key)")
end
if length(path) > 5
    println("  ...")
end

# Get path log-likelihood
ll = pathloglikelihood(sampler, t_max)
println("\nPath log-likelihood: $ll")
println()

# ===========================================================================
# Part 2: Replay function for computing log-likelihood of a given path
# ===========================================================================

"""
Compute the log-likelihood of a given path under the density-dependent birth-death model.
This "replays" the path through the simulation to compute the log-likelihood.
"""
function bd_path_logpdf(
    path::EventPath,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::Float64
    # Sanity checks
    if any(ev -> ev.time < 0 || ev.time > t_max, path)
        return -Inf
    end
    for i in 2:length(path)
        if path[i].time <= path[i-1].time
            return -Inf
        end
    end

    # Build sampler for replay (RNG is not used for sampling)
    rng = Xoshiro(0)
    sampler = SamplingContext(ClockKey, Float64, rng; path_likelihood=true)
    state = BDState(N0)

    # Density-dependent birth rate (must match simulate_bd)
    function birth_rate(N::Int)
        N <= 0 && return 0.0
        return λ_birth * N * max(0.0, 1.0 - N / K)
    end

    # Enable initial clocks
    rate = birth_rate(length(state.population))
    if rate > 0
        enable!(sampler, (:birth, 0), Exponential(inv(rate)))
    end
    for i in state.population
        enable!(sampler, (:death, i), Gamma(death_shape, death_scale))
    end

    # Replay each event
    for ev in path
        if !isenabled(sampler, ev.key)
            return -Inf  # Invalid path: firing disabled clock
        end

        fire!(sampler, ev.key, ev.time)

        # Handle event
        if ev.key[1] == :birth
            new_id = state.next_id
            state.next_id += 1
            push!(state.population, new_id)
            enable!(sampler, (:death, new_id), Gamma(death_shape, death_scale))
        elseif ev.key[1] == :death
            delete!(state.population, ev.key[2])
        end

        # Refresh birth clock with density-dependent rate
        rate = birth_rate(length(state.population))
        if rate > 0
            enable!(sampler, (:birth, 0), Exponential(inv(rate)))
        end
    end

    return pathloglikelihood(sampler, t_max)
end

# ===========================================================================
# Test Part 2: Verify replay gives same log-likelihood
# ===========================================================================

println("=" ^ 60)
println("Part 2: Replay function")
println("=" ^ 60)

ll_replay = bd_path_logpdf(path, t_max, λ_birth, K, death_shape, death_scale, N0)
println("Original log-likelihood: $ll")
println("Replayed log-likelihood: $ll_replay")
println("Match: $(isapprox(ll, ll_replay))")
println()

# ===========================================================================
# Part 3: Wrap as a Gen distribution over event paths
# ===========================================================================

println("Loading Gen...")
using Gen
println("Gen loaded.")

"""
Custom Gen distribution over birth-death event paths.
The value is an `EventPath` (Vector{BDEvent}).
"""
struct BDPathDist <: Gen.Distribution{EventPath} end
const bd_path_dist = BDPathDist()

# Sample a path from the BD process
function Gen.random(
    ::BDPathDist,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::EventPath
    rng = Random.default_rng()
    path, _ = simulate_bd(rng, t_max, λ_birth, K, death_shape, death_scale, N0)
    return path
end

# Log density of a given path
function Gen.logpdf(
    ::BDPathDist,
    path::EventPath,
    t_max::Float64,
    λ_birth::Float64,
    K::Float64,
    death_shape::Float64,
    death_scale::Float64,
    N0::Int,
)::Float64
    return bd_path_logpdf(path, t_max, λ_birth, K, death_shape, death_scale, N0)
end

# This is a continuous distribution over paths
Gen.is_discrete(::BDPathDist) = false

# No gradients provided (can be upgraded later)
Gen.has_output_grad(::BDPathDist) = false
Gen.has_argument_grads(::BDPathDist) = (false, false, false, false, false, false)

# ===========================================================================
# Test Part 3: Gen distribution standalone
# ===========================================================================

println("=" ^ 60)
println("Part 3: Gen distribution")
println("=" ^ 60)

# Sample from distribution
gen_path = Gen.random(bd_path_dist, t_max, λ_birth, K, death_shape, death_scale, N0)
println("Sampled path with $(length(gen_path)) events")

# Compute logpdf
gen_logp = Gen.logpdf(bd_path_dist, gen_path, t_max, λ_birth, K, death_shape, death_scale, N0)
println("log p(path | θ) = $gen_logp")
println()

# ===========================================================================
# Part 4: Using the distribution inside a Gen model
# ===========================================================================

"""
A Gen model that:
- Takes carrying capacity K as a fixed parameter
- Puts priors on birth rate and mean death time
- Samples a full event path as a single random choice

The density-dependent birth rate ensures stable population dynamics.
"""
@gen function bd_model(t_max::Float64, K::Float64, N0::Int)
    # Priors on parameters
    λ_birth = @trace(gamma(2.0, 1.0), :λ_birth)        # Gamma(shape=2, scale=1), mean=2
    mean_death = @trace(gamma(2.0, 2.0), :mean_death)  # Gamma(shape=2, scale=2), mean=4

    death_shape = 2.0
    death_scale = mean_death / death_shape

    # One random choice: the entire event path
    path = @trace(
        bd_path_dist(t_max, λ_birth, K, death_shape, death_scale, N0),
        :path,
    )

    return (λ_birth=λ_birth, mean_death=mean_death, path=path)
end

# ===========================================================================
# Test Part 4: Gen model
# ===========================================================================

println("=" ^ 60)
println("Part 4: Gen model")
println("=" ^ 60)

# Forward simulation
trace = Gen.simulate(bd_model, (t_max, K, N0))
retval = Gen.get_retval(trace)

println("Sampled λ_birth = $(retval.λ_birth)")
println("Sampled mean_death = $(retval.mean_death)")
println("Number of events = $(length(retval.path))")
println("Trace log-probability = $(Gen.get_score(trace))")
println()

# Conditioning on an observed path
println("Conditioning on observed path...")
obs = Gen.choicemap((:path, path))  # Use path from Part 1
trace_cond, logw = Gen.generate(bd_model, (t_max, K, N0), obs)

println("Conditioned λ_birth = $(trace_cond[:λ_birth])")
println("Conditioned mean_death = $(trace_cond[:mean_death])")
println("Importance weight = $logw")
println()

println("=" ^ 60)
println("Example complete!")
println("=" ^ 60)
