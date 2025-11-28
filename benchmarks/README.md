# CompetingClocks.jl Benchmarks

This directory contains a comprehensive benchmarking system for CompetingClocks.jl samplers. The benchmarks measure performance across different sampler types, workload scales, distribution types, and key strategies to help users choose the optimal sampler for their specific use case.

## Quick Start

### Running Benchmarks

**Production run with full output:**
```bash
cd benchmarks
julia --project=. src/run_all.jl
```

**Simple run without extras:**
```bash
cd benchmarks
julia --project=. src/benchmark_runner.jl
```

**Results:** Both scripts save benchmark results to `data/observations.csv`.

### Example Output

```
======================================================================
               CompetingClocks.jl Sampler Benchmarks
======================================================================

Configuration:
  Samplers: 3
  Conditions: 144
  Total combinations: 432
  Valid benchmarks: 336

[1/432] FirstToFire: n=10, churn=1, dist=exponential, keys=dense... ✓ 584ns
[2/432] FirstToFire: n=10, churn=1, dist=exponential, keys=sparse... ✓ 1.4μs
...

Summary:
  Completed: 336 benchmarks
  Skipped: 96 (incompatible combinations)
  Elapsed time: 1847.3s
  Output file: data/observations.csv
```

## Customizing Benchmarks

All benchmark configuration is centralized in `src/conditions.jl`. Edit the constants to customize which benchmarks run:

### Configuration Constants

```julia
# Number of enabled clocks to test
const N_ENABLED = [10, 100, 1_000, 10_000]

# Number of clocks to enable/disable per simulation step (churn rate)
const N_CHANGES = [1, 10, 100]

# Distribution types to test
const DISTRIBUTIONS = [:exponential, :gamma, :weibull]

# Key strategies: dense (1:n) vs sparse (random large integers)
const KEY_STRATEGIES = [:dense, :sparse]

# Sampler types to benchmark
const SAMPLERS = [
    FirstToFire{Int,Float64},
    DirectCall{Int,Float64},
    CombinedNextReaction{Int,Float64}
]
```

### Example: Quick Testing Configuration

For faster iteration during development:

```julia
const N_ENABLED = [10, 100]           # Just small scales
const N_CHANGES = [1, 5]              # Just low churn
const DISTRIBUTIONS = [:exponential]  # Just exponential
const KEY_STRATEGIES = [:dense]       # Just dense keys
const SAMPLERS = [FirstToFire{Int,Float64}]  # Just one sampler
```

This reduces the benchmark count from 432 to just 4 combinations.

### Example: Production Configuration

For comprehensive performance characterization:

```julia
const N_ENABLED = [10, 100, 1_000, 10_000, 100_000]
const N_CHANGES = [1, 10, 100, 1_000]
const DISTRIBUTIONS = [:exponential, :gamma, :weibull]
const KEY_STRATEGIES = [:dense, :sparse]
const SAMPLERS = [
    FirstToFire{Int,Float64},
    DirectCall{Int,Float64},
    CombinedNextReaction{Int,Float64}
]
```

## Understanding Results

### CSV Output Schema

The `data/observations.csv` file contains these columns:

| Column | Type | Description |
|--------|------|-------------|
| `sampler_type` | String | Sampler class name (FirstToFire, DirectCall, CombinedNextReaction) |
| `n_enabled` | Int | Number of enabled clocks in the sampler |
| `n_changes` | Int | Number of clocks enabled/disabled per step (churn rate) |
| `distributions` | String | Distribution type (exponential, gamma, weibull) |
| `key_strategy` | String | Key allocation strategy (dense, sparse) |
| `time_ns` | Int | Median execution time in nanoseconds |
| `memory_bytes` | Int | Median memory allocation in bytes |

### Interpreting Performance

**Time (time_ns):**
- Lower is better
- Measured in nanoseconds for single simulation step
- Includes: next(), fire!(), disable!(), and enable!() operations

**Memory (memory_bytes):**
- Lower is better
- Shows allocation churn per step, not peak memory
- High memory indicates frequent allocations

### Sampler Characteristics

**FirstToFire:**
- General-purpose sampler using binary heap
- Works with all distribution types
- Good performance across workload scales
- Recommended as default choice

**DirectCall:**
- Specialized for exponential distributions only
- Direct method without heap operations
- Can be faster for certain workloads
- Higher memory usage due to recalculation approach

**CombinedNextReaction:**
- Hybrid approach combining multiple techniques
- Works with all distribution types
- Often competitive with FirstToFire
- May excel at specific workload patterns

### Key Strategies

**Dense keys (`:dense`):**
- Uses sequential integers (1, 2, 3, ...)
- Best cache locality
- Reuses keys when re-enabling
- Represents simulations with fixed set of events

**Sparse keys (`:sparse`):**
- Uses random large integers
- Tests hash table performance
- Generates brand new keys each time
- Represents dynamic event creation scenarios

## Understanding Benchmark Conditions

### What is n_enabled?

The number of clocks that are enabled in the sampler at steady state. Think of this as the "pool size" of potential events.

**Example:** With `n_enabled=100`, the sampler maintains 100 active clocks that could fire.

### What is n_changes?

The "churn rate" - how many clocks get disabled and re-enabled per simulation step.

**Example:** With `n_changes=10`:
1. One clock fires (chosen by `next()` and disabled via `fire!()`)
2. Nine additional clocks are disabled
3. All ten clocks are re-enabled with new distributions

Higher churn rates stress the enable/disable operations more than the next() operation.

### Why Test Different Combinations?

Different simulation workloads have different characteristics:

- **Small n_enabled, low churn:** Simple discrete event simulations
- **Large n_enabled, low churn:** Many events, infrequent updates (e.g., epidemic models)
- **Large n_enabled, high churn:** Many events, frequent updates (e.g., chemical reactions)
- **Sparse keys:** Dynamically created events (e.g., agent-based models with births/deaths)

## Sampler Compatibility

Not all samplers work with all distribution types:

| Sampler | Exponential | Gamma | Weibull |
|---------|-------------|-------|---------|
| FirstToFire | ✓ | ✓ | ✓ |
| DirectCall | ✓ | ✗ | ✗ |
| CombinedNextReaction | ✓ | ✓ | ✓ |

**DirectCall** only supports exponential distributions because it uses the Direct Method algorithm, which requires constant hazard rates.

Incompatible combinations are automatically skipped during benchmarking.

## Module Structure

```
benchmarks/
├── src/
│   ├── Benchmarks.jl           # Main module definition
│   ├── conditions.jl           # Configuration constants and types
│   ├── measure.jl              # Benchmarking functions
│   ├── benchmark_runner.jl     # Simple benchmark runner
│   ├── run_all.jl             # Production runner with full UX
│   └── minimal_benchmark.jl   # Original prototype (superseded)
├── data/
│   └── observations.csv       # Benchmark results
├── Project.toml               # Julia package configuration
└── README.md                  # This file
```

### Module Components

**`Benchmarks.jl`**
- Module entry point
- Exports key functions and constants

**`conditions.jl`**
- `BenchmarkCondition` struct
- Configuration constants (edit these to customize)
- `generate_conditions()` - Creates all valid combinations
- `is_compatible()` - Checks sampler/condition compatibility

**`measure.jl`**
- `create_distribution()` - Instantiates distribution objects
- `setup_sampler()` - Initializes sampler with condition
- `benchmark_step!()` - Performs one simulation step
- `benchmark_config()` - Runs complete benchmark with BenchmarkTools

## Advanced Usage

### Using the Module in Your Code

```julia
include("src/Benchmarks.jl")
using .Benchmarks

# Generate conditions
conditions = generate_conditions()

# Filter to specific conditions
small_conditions = filter(c -> c.n_enabled <= 100, conditions)

# Run single benchmark
sampler_type = FirstToFire{Int,Float64}
cond = BenchmarkCondition(100, 10, :exponential, :dense)
time_ns, mem_bytes = benchmark_config(sampler_type, cond)

println("Time: $(time_ns)ns, Memory: $(mem_bytes) bytes")
```

### Creating Custom Benchmarks

```julia
include("src/Benchmarks.jl")
using .Benchmarks
using CSV, DataFrames

# Define custom conditions
custom_conditions = [
    BenchmarkCondition(50, 5, :exponential, :dense),
    BenchmarkCondition(200, 20, :gamma, :sparse),
]

# Run benchmarks
results = DataFrame(...)
for cond in custom_conditions
    time_ns, mem = benchmark_config(FirstToFire{Int,Float64}, cond)
    push!(results, (...))
end

CSV.write("custom_results.csv", results)
```

### Modifying Benchmark Parameters

The `benchmark_config()` function uses `@benchmark` from BenchmarkTools with:
- `samples=100` - Number of timing samples

To adjust these parameters, edit `src/measure.jl`:

```julia
result = @benchmark benchmark_step!(
    $sampler, $enabled_keys, $(cond.key_strategy), $(cond.n_changes),
    $dist_type, $rng, $when
) samples=100 evals=1  # Adjust here
```

## Tips for Accurate Benchmarking

1. **Close other applications** to minimize system noise
2. **Disable CPU frequency scaling** for consistent results
3. **Run multiple times** and compare results for variability
4. **Use production configuration** to identify trends across scales
5. **Watch for thermal throttling** on long benchmark runs

## Analyzing Results

### Loading Results

```julia
using CSV, DataFrames

results = CSV.read("data/observations.csv", DataFrame)

# Filter to specific sampler
ftf_results = filter(row -> row.sampler_type == "FirstToFire", results)

# Compare samplers at specific scale
scale_100 = filter(row -> row.n_enabled == 100, results)

# Find fastest sampler per condition
using DataFrames
fastest = combine(groupby(results, [:n_enabled, :n_changes, :distributions, :key_strategy]),
                  :time_ns => minimum => :min_time,
                  :sampler_type => (x -> x[argmin(results.time_ns)]) => :fastest_sampler)
```

### Visualization Ideas

- **Scaling plots:** Time vs n_enabled for each sampler
- **Heatmaps:** Performance across n_enabled × n_changes
- **Comparisons:** Sampler performance ratios
- **Memory trends:** Allocation patterns across scales

## Troubleshooting

### Benchmarks Taking Too Long

**Reduce the parameter space:**
```julia
const N_ENABLED = [10, 100]  # Remove large scales
const N_CHANGES = [1, 10]    # Remove high churn
```

### Out of Memory

**Test smaller scales first:**
```julia
const N_ENABLED = [10, 100, 1_000]  # Skip 10_000
```

### Inconsistent Results

**Increase sample count in `measure.jl`:**
```julia
) samples=200  # Up from 100
```

## Contributing

When adding new samplers or benchmark conditions:

1. Add sampler type to `SAMPLERS` in `conditions.jl`
2. Add compatibility check to `is_compatible()` if needed
3. Test with small configuration first
4. Document any special requirements

## References

- [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) - Julia benchmarking framework
- [CompetingClocks.jl Documentation](../docs) - Main package documentation
- [Benchmark Requirements](../docs/dev/benchmark_rqmts.md) - Original requirements document
- [Implementation Plan](../docs/dev/benchmark_plan.md) - Detailed implementation steps
