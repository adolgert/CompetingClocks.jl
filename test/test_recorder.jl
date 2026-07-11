using SafeTestsets
using Test
using CompetingClocks
using CompetingClocks: with_recorder, close_record!, recorded_firings
using Distributions
using Random: Xoshiro


# A five-machine failure/repair model driven through a SamplingContext. Each
# machine has a failure clock while working and a repair clock while broken;
# firing one enables the other. Keys are (:fail, i) / (:repair, i). Returns the
# ordered list of (key, when) the driver itself fired — an independent record of
# what happened, captured OUTSIDE the recorder under test.
function run_machine_repair!(ctx, nmachines, horizon, fail_dist, repair_dist)
    for i in 1:nmachines
        enable!(ctx, (:fail, i), fail_dist(i))
    end
    fired = Tuple{Tuple{Symbol,Int},Float64}[]
    while true
        when, which = next(ctx)
        when > horizon && break
        fire!(ctx, which, when)
        push!(fired, (which, when))
        kind, i = which
        if kind === :fail
            enable!(ctx, (:repair, i), repair_dist(i))
        else
            enable!(ctx, (:fail, i), fail_dist(i))
        end
    end
    return fired
end


# Plain @testset (not @safetestset) so the file-scope helper above is in scope,
# matching the mixed style of test_context.jl.
@testset "recorder: every recorded firing satisfies when == te + invlogccdf(dist, log u) to 1e-9 on both FirstReaction and CombinedNextReaction" begin
    # Weibull failures (non-exponential, so age genuinely matters) racing
    # exponential repairs. The identity is a contract obligation the recorder
    # must honor no matter which sampler produced the firing time.
    fail_dist(i) = Weibull(1.7, 1.0 + 0.1 * i)
    repair_dist(i) = Exponential(0.5)
    horizon = 40.0

    for method in (FirstReactionMethod(), NextReactionMethod())
        ctx = SamplingContext(SamplerBuilder(Tuple{Symbol,Int}, Float64; method=method),
                              Xoshiro(20260709))
        ctx, rec = with_recorder(ctx)
        run_machine_repair!(ctx, 5, horizon, fail_dist, repair_dist)
        close_record!(rec, horizon)

        firings = recorded_firings(rec)
        # A run this long must actually exercise many firings, else the identity
        # would be asserted on an empty set and prove nothing.
        @test length(firings) > 20
        for fr in firings
            @test isfinite(fr.u) && 0.0 <= fr.u <= 1.0
            # The identity, in the log-survival space where it lives.
            @test isapprox(fr.when, fr.te + invlogccdf(fr.distribution, fr.logu); atol=1e-9)
            # And equivalently through u itself.
            @test isapprox(fr.when, fr.te + invlogccdf(fr.distribution, log(fr.u)); atol=1e-9)
        end
    end
end


@testset "recorder: the recorded firing sequence of keys and times equals the driver's own independently captured firing history" begin
    fail_dist(i) = Weibull(1.4, 2.0)
    repair_dist(i) = Exponential(0.7)
    horizon = 30.0

    for method in (FirstReactionMethod(), NextReactionMethod())
        ctx = SamplingContext(SamplerBuilder(Tuple{Symbol,Int}, Float64; method=method),
                              Xoshiro(555))
        ctx, rec = with_recorder(ctx)
        # The driver returns exactly the (key, when) pairs it fired — an oracle
        # built without consulting the recorder.
        driver_history = run_machine_repair!(ctx, 5, horizon, fail_dist, repair_dist)
        close_record!(rec, horizon)

        firings = recorded_firings(rec)
        @test length(firings) == length(driver_history)
        for (fr, (key, when)) in zip(firings, driver_history)
            @test fr.clock == key
            @test fr.when == when
        end
    end
end


@safetestset recorder_left_shifted_enabling =
    "recorder: a left-shifted enabling (te < when) records the total-lifetime uniform u == ccdf(dist, when - te), not one conditioned on survival to the enable call" begin
    using CompetingClocks
    using CompetingClocks: with_recorder, close_record!, recorded_firings
    using Distributions: Weibull, ccdf, invlogccdf
    using Random: Xoshiro

    # One clock, enabled at simulation time 0 but with its distribution's zero
    # shifted 1.5 into the past (te = -1.5): the clock is already aged 1.5 when
    # observed. The recorded u must be the survival of the TOTAL lifetime from
    # te, so u == ccdf(dist, when - te) with when - te > 1.5.
    dist = Weibull(1.7, 2.0)
    age = 1.5

    for method in (FirstReactionMethod(), NextReactionMethod())
        ctx = SamplingContext(SamplerBuilder(Int, Float64; method=method), Xoshiro(31415))
        ctx, rec = with_recorder(ctx)
        enable!(ctx, 1, dist, -age)     # relative_te = -1.5  ⇒  te = -1.5 < when = 0
        when, which = next(ctx)
        fire!(ctx, which, when)
        close_record!(rec, when + 10.0)

        fr = recorded_firings(rec)[1]
        @test fr.te == -age
        @test fr.when - fr.te > age                # aged past the shift
        # The stored uniform is the total-lifetime survival from te, NOT the
        # survival conditioned on having reached the enable time (which would be
        # ccdf(dist, when - te) / ccdf(dist, age)).
        @test isapprox(fr.u, ccdf(dist, fr.when - fr.te); atol=1e-12)
        @test isapprox(fr.when, fr.te + invlogccdf(dist, fr.logu); atol=1e-9)
        # Guard: the conditioned-uniform mistake would give a strictly larger
        # number, so the test can actually see the difference.
        conditioned = ccdf(dist, fr.when - fr.te) / ccdf(dist, age)
        @test conditioned > fr.u
    end
end


@safetestset recorder_reset_and_close =
    "recorder: reset! clears firings and horizon so a reused recorder starts empty, and close_record! stamps the horizon" begin
    using CompetingClocks
    using CompetingClocks: TrajectoryRecorder, close_record!, recorded_firings, isclosed, horizon, reset!
    using Distributions: Exponential
    using Random: Xoshiro

    rng = Xoshiro(9001)
    rec = TrajectoryRecorder{Int,Float64}()
    @test !isclosed(rec)

    enable!(rec, 1, Exponential(1.0), 0.0, 0.0)
    fire!(rec, 1, 0.8)
    close_record!(rec, 5.0)
    @test length(recorded_firings(rec)) == 1
    @test isclosed(rec)
    @test horizon(rec) == 5.0

    reset!(rec)
    @test isempty(recorded_firings(rec))
    @test !isclosed(rec)
    @test horizon(rec) == 0.0
    # Reusable after reset.
    enable!(rec, 2, Exponential(2.0), 0.0, 0.0)
    fire!(rec, 2, 1.1)
    @test length(recorded_firings(rec)) == 1
    @test recorded_firings(rec)[1].clock == 2
end
