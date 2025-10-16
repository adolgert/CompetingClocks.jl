using BenchmarkTools
using CompetingClocks
using CSV, DataFrames
using Distributions
using Random

# Import from Base.Iterators
import Base.Iterators: take, flatten

# Step 3: Create distribution based on type
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

# Step 3 & 4: Create sampler with parameters for distribution and key strategy
function setup_minimal_sampler(n_enabled::Int, dist_type::Symbol, key_strategy::Symbol)
    sampler = FirstToFire{Int,Float64}()
    rng = MersenneTwister(42)

    # Generate keys based on strategy
    keys = if key_strategy == :dense
        collect(1:n_enabled)
    elseif key_strategy == :sparse
        rand(rng, 1:typemax(Int32), n_enabled)
    else
        error("Unknown key strategy: $key_strategy")
    end

    # Enable all clocks
    for key in keys
        dist = create_distribution(dist_type, rng)
        enable!(sampler, key, dist, 0.0, 0.0, rng)
    end

    return sampler, keys, rng, dist_type
end

# Step 2: One complete simulation step with churn (multiple enable/disable)
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

    if key_strategy == :dense
        all_keys_to_enable = flatten([[what_fire], disable_keys])
    elseif key_strategy == :sparse
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

    # 4. Re-enable all n_changes clocks with new distributions
    for k in first(all_keys_to_enable, n_changes)
        dist = create_distribution(dist_type, rng)
        enable!(sampler, k, dist, when_fire, when_fire, rng)
        push!(enabled_keys, k)
    end

    return when_fire
end

# Run the benchmark
function main()
    # Steps 2, 3, 4: Test all combinations
    results = DataFrame(
        sampler_type = String[],
        n_enabled = Int[],
        n_changes = Int[],
        distributions = String[],
        key_strategy = String[],
        time_ns = Int[],
        memory_bytes = Int[]
    )

    n_enabled = 10
    n_changes_values = [1, 5]
    dist_types = [:exponential, :gamma, :weibull]
    key_strategies = [:dense, :sparse]

    total = length(n_changes_values) * length(dist_types) * length(key_strategies)
    done = 0

    for n_changes in n_changes_values
        for dist_type in dist_types
            for key_strategy in key_strategies
                done += 1
                println("[$done/$total] Running benchmark: n_changes=$n_changes, dist=$dist_type, keys=$key_strategy")

                sampler, keys, rng, dist_type_val = setup_minimal_sampler(n_enabled, dist_type, key_strategy)
                when = 0.0
                keys = Set(keys)

                result = @benchmark benchmark_step!(
                    $sampler, $keys, $key_strategy, $n_changes, $dist_type_val, $rng, $when
                ) samples=100

                push!(results, (
                    "FirstToFire",
                    n_enabled,
                    n_changes,
                    string(dist_type),
                    string(key_strategy),
                    Int(round(median(result.times))),
                    Int(round(median(result.memory)))
                ))
            end
        end
    end

    # Save to CSV
    mkpath("data")
    CSV.write("data/observations.csv", results)

    println("\nResults:")
    println(results)
    println("\nSaved to data/observations.csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end