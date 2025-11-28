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

@safetestset keyedkeep_copy = "KeyedKeepPrefixSearch copy!" begin
    using CompetingClocks

    prefix_src = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    src = CompetingClocks.KeyedKeepPrefixSearch{String,typeof(prefix_src)}(prefix_src)
    src["a"] = 1.0
    src["b"] = 2.0
    src["c"] = 3.0

    prefix_dst = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    dst = CompetingClocks.KeyedKeepPrefixSearch{String,typeof(prefix_dst)}(prefix_dst)
    dst["x"] = 0.0
    dst["y"] = 0.0
    dst["z"] = 0.0

    copy!(dst, src)

    @test dst["a"] == 1.0
    @test dst["b"] == 2.0
    @test dst["c"] == 3.0
    @test length(dst) == 3
end

@safetestset keyedkeep_key_type = "KeyedKeepPrefixSearch key_type" begin
    using CompetingClocks
    using CompetingClocks: key_type

    prefix = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    kp_string = CompetingClocks.KeyedKeepPrefixSearch{String,typeof(prefix)}(prefix)
    @test key_type(kp_string) == String

    prefix2 = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    kp_int = CompetingClocks.KeyedKeepPrefixSearch{Int,typeof(prefix2)}(prefix2)
    @test key_type(kp_int) == Int
end

@safetestset keyedkeep_prefixenabled = "PrefixEnabled in and eltype" begin
    using CompetingClocks
    using CompetingClocks: enabled, isenabled

    prefix = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    kp = CompetingClocks.KeyedKeepPrefixSearch{String,typeof(prefix)}(prefix)
    kp["a"] = 1.0
    kp["b"] = 0.0  # disabled (zero)
    kp["c"] = 3.0

    en = enabled(kp)

    # Test `in` for PrefixEnabled
    @test "a" in en
    @test !("b" in en)  # zero value means not enabled
    @test "c" in en
    @test !("d" in en)  # not present

    # Test eltype
    @test eltype(typeof(en)) == String
end

@safetestset keyedremoval_key_time_type = "KeyedRemovalPrefixSearch key_type and time_type" begin
    using CompetingClocks
    using CompetingClocks: key_type, time_type

    prefix = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    kp = CompetingClocks.KeyedRemovalPrefixSearch{String,typeof(prefix)}(prefix)

    @test key_type(kp) == String
    @test time_type(kp) == Float64
end

@safetestset keyedremoval_update_existing = "KeyedRemovalPrefixSearch update existing clock" begin
    using CompetingClocks

    prefix = CompetingClocks.BinaryTreePrefixSearch{Float64}()
    kp = CompetingClocks.KeyedRemovalPrefixSearch{String,typeof(prefix)}(prefix)

    # Add initial values
    kp["a"] = 1.0
    kp["b"] = 2.0

    @test kp["a"] == 1.0
    @test kp["b"] == 2.0

    # Update existing clock (not delete and re-add, just update)
    kp["a"] = 5.0
    kp["b"] = 6.0

    @test kp["a"] == 5.0
    @test kp["b"] == 6.0
    @test length(kp) == 2
end
