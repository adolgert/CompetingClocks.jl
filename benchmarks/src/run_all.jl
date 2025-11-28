include("Benchmarks.jl")
using .Benchmarks
using CSV, DataFrames

"""
Run the complete benchmark suite for CompetingClocks.jl samplers.

This script benchmarks all sampler types across all valid condition combinations,
providing comprehensive performance data for sampler selection.
"""
function main()
    println("=" ^ 70)
    println(" " ^ 15, "CompetingClocks.jl Sampler Benchmarks")
    println("=" ^ 70)
    println()

    conditions = generate_conditions()
    samplers = construct_samplers()

    # Calculate how many benchmarks will actually run
    total_combinations = length(samplers) * length(conditions)
    actual_benchmarks = sum(is_compatible(typeof(s), c) for s in samplers, c in conditions)

    println("Configuration:")
    println("  Samplers: $(length(samplers))")
    println("  Conditions: $(length(conditions))")
    println("  Total combinations: $total_combinations")
    println("  Valid benchmarks: $actual_benchmarks")
    println()
    println("Sampler types:")
    for s in samplers
        println("  • $(sampler_name(typeof(s)))")
    end
    println()
    println("Configuration ranges:")
    println("  N_ENABLED: $(N_ENABLED)")
    println("  N_CHANGES: $(N_CHANGES)")
    println("  DISTRIBUTIONS: $(DISTRIBUTIONS)")
    println("  KEY_STRATEGIES: $(KEY_STRATEGIES)")
    println()
    println("=" ^ 70)
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

    done = 0
    skipped = 0
    start_time = time()

    for sampler in samplers
        sampler_type = typeof(sampler)
        for cond in conditions
            done += 1

            # Skip incompatible sampler+condition combinations
            if !is_compatible(sampler_type, cond)
                skipped += 1
                println("[$done/$total_combinations] SKIP: $(nameof(sampler_type)) + $(cond.distributions) (incompatible)")
                continue
            end

            print("[$done/$total_combinations] $(sampler_name(sampler_type)): ")
            print("n=$(cond.n_enabled), churn=$(cond.n_changes), ")
            print("dist=$(cond.distributions), keys=$(cond.key_strategy)... ")
            flush(stdout)

            try
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

                # Format time nicely
                if time_ns < 1_000
                    time_str = "$(time_ns)ns"
                elseif time_ns < 1_000_000
                    time_str = "$(round(time_ns / 1_000, digits=1))μs"
                elseif time_ns < 1_000_000_000
                    time_str = "$(round(time_ns / 1_000_000, digits=1))ms"
                else
                    time_str = "$(round(time_ns / 1_000_000_000, digits=2))s"
                end

                println("✓ $time_str")
            catch e
                println("✗ ERROR")
                println("  Error: $e")
                continue
            end
        end
    end

    elapsed = time() - start_time

    # Save results
    mkpath("data")
    output_file = "data/observations.csv"
    CSV.write(output_file, results)

    println()
    println("=" ^ 70)
    println(" " ^ 25, "Benchmarking Complete")
    println("=" ^ 70)
    println()
    println("Summary:")
    println("  Completed: $(nrow(results)) benchmarks")
    println("  Skipped: $skipped (incompatible combinations)")
    println("  Elapsed time: $(round(elapsed, digits=1))s")
    println("  Output file: $output_file")
    println()

    # Show sample results
    if nrow(results) > 0
        println("Sample results (first 10):")
        println(first(results, min(10, nrow(results))))
        println()

        # Show fastest and slowest for each sampler at largest scale
        println("Performance highlights:")
        for sampler_name in unique(results.sampler_type)
            sampler_results = filter(row -> row.sampler_type == sampler_name, results)
            if nrow(sampler_results) > 0
                fastest = sampler_results[argmin(sampler_results.time_ns), :]
                slowest = sampler_results[argmax(sampler_results.time_ns), :]

                println("  $sampler_name:")
                println("    Fastest: $(fastest.time_ns)ns (n=$(fastest.n_enabled), churn=$(fastest.n_changes), $(fastest.distributions), $(fastest.key_strategy))")
                println("    Slowest: $(slowest.time_ns)ns (n=$(slowest.n_enabled), churn=$(slowest.n_changes), $(slowest.distributions), $(slowest.key_strategy))")
            end
        end
    end

    println()
    println("=" ^ 70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
