using CompetingClocks
using CompetingClocks: DirectCallExplicit, KeyedRemovalPrefixSearch, KeyedKeepPrefixSearch,
                        BinaryTreePrefixSearch, CumSumPrefixSearch

"""
Represents one benchmark condition configuration.
"""
struct BenchmarkCondition
    n_enabled::Int
    n_changes::Int
    distributions::Symbol
    key_strategy::Symbol
end

# User-editable configuration constants
const N_ENABLED = [10, 100, 1_000, 10_000]
const N_CHANGES = [1, 10, 100]
const DISTRIBUTIONS = [:exponential, :gamma, :weibull]
const KEY_STRATEGIES = [:dense, :sparse]
function construct_samplers()
    K = Int
    T = Float64
    return [
        FirstToFire{K,T}(),
        DirectCallExplicit(K, T, KeyedRemovalPrefixSearch, BinaryTreePrefixSearch),
        DirectCallExplicit(K, T, KeyedKeepPrefixSearch, BinaryTreePrefixSearch),
        DirectCallExplicit(K, T, KeyedRemovalPrefixSearch, CumSumPrefixSearch),
        DirectCallExplicit(K, T, KeyedKeepPrefixSearch, CumSumPrefixSearch),
        FirstReaction{K,T}(),
        CombinedNextReaction{K,T}()
    ]
end

"""
Generate all valid benchmark conditions.

Filters out invalid combinations where n_changes > n_enabled.
"""
function generate_conditions()
    conditions = BenchmarkCondition[]

    for n in N_ENABLED, m in N_CHANGES, d in DISTRIBUTIONS, k in KEY_STRATEGIES
        # Skip invalid conditions where we try to change more clocks than exist
        if m < n
            push!(conditions, BenchmarkCondition(n, m, d, k))
        end
        push!(conditions, BenchmarkCondition(n, n, d, k))
    end

    return conditions
end

"""
Check if a sampler type is compatible with a given condition.

DirectCall only works with exponential distributions.
"""
function is_compatible(sampler_type::Type, cond::BenchmarkCondition)
    # DirectCall can only handle exponential distributions
    if nameof(sampler_type) == :DirectCall && cond.distributions != :exponential
        return false
    end

    return true
end
