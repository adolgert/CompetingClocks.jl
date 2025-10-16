include("Benchmarks.jl")
using .Benchmarks
using CSV, DataFrames

"""
Run benchmarks for all conditions and save results to CSV.
"""
function main()
    conditions = generate_conditions()
    samplers = construct_samplers()

    println("Running $(length(samplers)) sampler(s) × $(length(conditions)) conditions = $(length(samplers) * length(conditions)) benchmarks")
    println()

    # Create results DataFrame
    results = DataFrame(
        sampler_type = String[],
        n_enabled = Int[],
        n_changes = Int[],
        distributions = String[],
        key_strategy = String[],
        time_ns = Int[],
        memory_bytes = Int[]
    )

    total = length(samplers) * length(conditions)
    done = 0

    for sampler in samplers
        sampler_type = typeof(sampler)
        for cond in conditions
            done += 1

            # Skip incompatible sampler+condition combinations
            if !is_compatible(sampler_type, cond)
                println("[$done/$total] SKIP: $(nameof(sampler_type)) (incompatible with $(cond.distributions))")
                continue
            end

            println("[$done/$total] Running: $(sampler_name(sampler_type)), n_enabled=$(cond.n_enabled), n_changes=$(cond.n_changes), dist=$(cond.distributions), keys=$(cond.key_strategy)")

            time_ns, mem_bytes = benchmark_config(sampler, cond)

            push!(results, (
                sampler_name(sampler_type),
                cond.n_enabled,
                cond.n_changes,
                string(cond.distributions),
                string(cond.key_strategy),
                time_ns,
                mem_bytes
            ))
        end
    end

    # Save to CSV
    mkpath("data")
    CSV.write("data/observations.csv", results)

    println()
    println("Results:")
    println(results)
    println()
    println("Saved to data/observations.csv")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
