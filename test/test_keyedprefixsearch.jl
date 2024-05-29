using SafeTestsets

module KeyedPrefixScanFixture
using Random
using Distributions
using CompetingClocks
using Test

@enum ScanActions EnableNew EnableOld Disable Empty CheckAll

function scan_gymnasium(keyed_scan)
    rng = Xoshiro(4209422)
    enabled = Dict{Int64,Float64}()
    disabled = Set{Int64}()
    prob = Dict(EnableNew => 0.4, EnableOld => 0.2, Disable => 0.4, Empty => 0.02, CheckAll => 0.1)
    p = [prob[a] for a in instances(ScanActions)]
    p = p ./ sum(p)
    action_dist = Categorical(p)
    min_hazard = 1.9
    key_cnt = 0

    for trial_idx in 1:1000
        action = instances(ScanActions)[rand(rng, action_dist)]
        if action == EnableNew
            key_cnt += 1
            value = min_hazard + rand(rng)
            keyed_scan[string(key_cnt)] = value
            enabled[key_cnt] = value
        elseif action == EnableOld
            if length(disabled) > 0
                reenable = pop!(disabled)
                value = min_hazard + rand(rng)
                keyed_scan[string(reenable)] = value
                enabled[reenable] = value
            end
        elseif action == Disable
            if length(enabled) > 0
                dis, value = pop!(enabled)
                delete!(keyed_scan, string(dis))
                push!(disabled, dis)
            end
        elseif action == Empty
            empty!(keyed_scan)
            empty!(enabled)
            empty!(disabled)
        elseif action == CheckAll
            total = sum!(keyed_scan)
            if total > 0.0
                @test total > min_hazard * length(enabled)
                probed = Set{Int64}()
                # Walk slowly enough that we must see all keys.
                for probe in zero(total):(0.5 * min_hazard):total
                    key, hazard = CompetingClocks.choose(keyed_scan, probe)
                    push!(probed, parse(Int, key))
                end
                @test probed == Set(keys(enabled))
            else
                @test length(enabled) == 0
            end
        else
            @test action ∈ instances(Actions)
        end
    end
end
end


@safetestset keyedremoval_smoke = "KeyedRemoval smoke" begin
    using ..KeyedPrefixScanFixture: scan_gymnasium
    using Random
    using Distributions
    using CompetingClocks
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    keyed_scan = CompetingClocks.KeyedRemovalPrefixSearch{String,typeof(prefix_tree)}(prefix_tree)
    scan_gymnasium(keyed_scan)
end


@safetestset keyedkeep_smoke = "KeyedRemoval smoke" begin
    using ..KeyedPrefixScanFixture: scan_gymnasium
    using Random
    using Distributions
    using CompetingClocks
    prefix_tree = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    keyed_scan = CompetingClocks.KeyedKeepPrefixSearch{String,typeof(prefix_tree)}(prefix_tree)
    scan_gymnasium(keyed_scan)
end
