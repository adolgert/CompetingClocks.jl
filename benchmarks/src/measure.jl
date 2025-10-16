using BenchmarkTools
using CompetingClocks
using Distributions
using Random

# Import from Base.Iterators
import Base.Iterators: take, flatten

"""
Create a distribution instance based on the distribution type symbol.
"""
function create_distribution(dist_type::Symbol, rng::AbstractRNG)
    if dist_type == :exponential
        return Exponential(1.0)
    elseif dist_type == :gamma
        return Gamma(2.0, 1.0)
    elseif dist_type == :weibull
        return Weibull(2.0, 1.0)
    else
        error("Unknown distribution type: $dist_type")
    end
end

"""
Set up a sampler with the given condition parameters.

Returns (sampler, enabled_keys, rng, dist_type) where enabled_keys is a Set
that tracks which keys are currently enabled.
"""
function setup_sampler(sampler, cond::BenchmarkCondition)
    reset!(sampler)
    rng = MersenneTwister(42)

    # Generate keys based on strategy
    keys = if cond.key_strategy == :dense
        collect(1:cond.n_enabled)
    elseif cond.key_strategy == :sparse
        rand(rng, 1:typemax(Int32), cond.n_enabled)
    else
        error("Unknown key strategy: $(cond.key_strategy)")
    end

    # Enable all clocks
    for key in keys
        dist = create_distribution(cond.distributions, rng)
        enable!(sampler, key, dist, 0.0, 0.0, rng)
    end

    return Set(keys), rng, cond.distributions
end

"""
Perform one complete simulation step with churn (multiple enable/disable operations).

This function mutates both the sampler and the enabled_keys set to keep them synchronized.
"""
function benchmark_step!(sampler, enabled_keys, key_strategy, n_changes, dist_type, rng, when)
    # 1. Get next event
    (when_fire, what_fire) = next(sampler, when, rng)

    # 2. Fire it (disables the clock)
    fire!(sampler, what_fire, when_fire)
    delete!(enabled_keys, what_fire)

    # 3. Select n_changes-1 additional keys to disable
    disable_keys = collect(take(filter(k -> k != what_fire, enabled_keys), n_changes - 1))
    for k in disable_keys
        disable!(sampler, k, when_fire)
        delete!(enabled_keys, k)
    end

    # 4. Determine which keys to re-enable
    if key_strategy == :dense
        all_keys_to_enable = flatten([[what_fire], disable_keys])
    elseif key_strategy == :sparse
        # Generate brand new random keys for sparse strategy
        all_keys_to_enable = Int[]
        sizehint!(all_keys_to_enable, n_changes)
        while length(all_keys_to_enable) < n_changes
            batch = rand(rng, 1:typemax(Int32), 16)
            valid_keys = filter(k -> k ∉ enabled_keys && k ∉ all_keys_to_enable, batch)
            needed = n_changes - length(all_keys_to_enable)
            append!(all_keys_to_enable, first(valid_keys, min(needed, length(valid_keys))))
        end
    else
        error("Unknown key strategy: $key_strategy")
    end

    # 5. Re-enable all n_changes clocks with new distributions
    for k in first(all_keys_to_enable, n_changes)
        dist = create_distribution(dist_type, rng)
        enable!(sampler, k, dist, when_fire, when_fire, rng)
        push!(enabled_keys, k)
    end

    return when_fire
end

"""
Run a complete benchmark for a given sampler type and condition.

Returns (median_time_ns, median_memory_bytes).
"""
function benchmark_config(sampler, cond::BenchmarkCondition)
    # Set up the sampler
    enabled_keys, rng, dist_type = setup_sampler(sampler, cond)
    when = 0.0

    # Run the benchmark
    result = @benchmark benchmark_step!(
        $sampler, $enabled_keys, $(cond.key_strategy), $(cond.n_changes),
        $dist_type, $rng, $when
    ) samples=100

    return Int(round(median(result.times))), Int(round(median(result.memory)))
end
