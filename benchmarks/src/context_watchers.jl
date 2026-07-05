# Benchmark: SamplingContext watcher fan-out parity.
# Birth-death simulation, Int keys / Float64 time, ~100_000 events.
# Three configs: (a) no watchers, (b) recording=true, (c) path_likelihood=true.
# Reports min-of-5 wall time and allocations, plus enable! inference cleanliness.

using CompetingClocks
using CompetingClocks: SamplingContext, SamplerBuilder, enable!, disable!, fire!, next, reset!,
    FirstToFireMethod
using Random
using Distributions
using InteractiveUtils

# Deterministic birth-death. Two clocks: 1 = birth, 2 = death.
# Population-dependent rates; immigration keeps population from vanishing.
function run_bd!(ctx, nevents::Int)
    reset!(ctx)
    pop = 10
    b = 0.6; d = 0.5; immig = 1.0
    # initial enable
    enable!(ctx, 1, Exponential(1.0 / (b * pop + immig)))
    enable!(ctx, 2, Exponential(1.0 / (d * pop)))
    fired = 0
    while fired < nevents
        when, which = next(ctx)
        fire!(ctx, which, when)
        if which == 1
            pop += 1
        else
            pop = max(pop - 1, 1)
        end
        # `which` was auto-disabled by fire!; disable the survivor, then re-enable both.
        other = which == 1 ? 2 : 1
        disable!(ctx, other)
        enable!(ctx, 1, Exponential(1.0 / (b * pop + immig)))
        enable!(ctx, 2, Exponential(1.0 / (d * pop)))
        fired += 1
    end
    return pop
end

function make_ctx(config, seed)
    rng = Xoshiro(seed)
    if config == :none
        builder = SamplerBuilder(Int64, Float64; method=FirstToFireMethod())
    elseif config == :recording
        builder = SamplerBuilder(Int64, Float64; recording=true, method=FirstToFireMethod())
    elseif config == :path
        builder = SamplerBuilder(Int64, Float64; path_likelihood=true)
    else
        error("unknown config $config")
    end
    return SamplingContext(builder, rng)
end

function bench(config, nevents; runs=5)
    # warmup / compile
    ctx = make_ctx(config, 12345)
    run_bd!(ctx, 1000)
    tmin = Inf
    amin = typemax(Int)
    for _ in 1:runs
        ctx = make_ctx(config, 12345)
        a = @allocated run_bd!(ctx, nevents)
        ctx = make_ctx(config, 12345)
        t = @elapsed run_bd!(ctx, nevents)
        tmin = min(tmin, t)
        amin = min(amin, a)
    end
    return tmin, amin
end

function inference_report()
    ctx = make_ctx(:none, 999)
    enable!(ctx, 1, Exponential(0.5), 0.0)
    res = @code_typed enable!(ctx, 1, Exponential(0.5), 0.0)
    ci = res.first          # CodeInfo
    rt = res.second         # inferred return type
    branch_count = count(e -> isa(e, Core.GotoIfNot), ci.code)
    # scan for any non-concrete (Any) inferred ssa value
    any_any = any(t -> t === Any, ci.ssavaluetypes)
    return branch_count, rt, any_any, length(ci.code)
end

const NEV = 100_000

println("=== SamplingContext watcher benchmark ===")
println("julia ", VERSION)
for config in (:none, :recording, :path)
    t, a = bench(config, NEV)
    println(rpad(string(config), 12), "  time(min5)=", round(t * 1e3, digits=3), " ms",
        "   alloc(min5)=", a, " bytes")
end
bc, rt, anyany, ncode = inference_report()
println("--- enable! inference (config=none) ---")
println("GotoIfNot branches: ", bc)
println("return type: ", rt)
println("any Any in ssavaluetypes: ", anyany)
println("typed code stmts: ", ncode)
