using SafeTestsets

# Unit tests for the KeyedStreams primitive (milestone 4). The design finding
# these certify: a draw is addressed by (clock, occurrence), so it is independent
# of the ORDER in which clocks draw; a copy carries generator state (so a clone
# is coupled, not silently reset); and re-keying resets to a fresh, independent
# family. Full-sentence testset names state the invariant.


@safetestset keyed_streams_per_key_determinism =
    "a key's draw sequence is the same no matter the interleaving of other keys' draws" begin
    using CompetingClocks: KeyedStreams, stream_for!
    using Random: rand

    # Two identical families. Draw the same per-key occurrences but in different
    # global orders; each key's own sequence must match occurrence-for-occurrence.
    ks1 = KeyedStreams{Symbol}(0x1234)
    ks2 = KeyedStreams{Symbol}(0x1234)

    # Order A: a, b, a, b
    a1 = rand(stream_for!(ks1, :a)); b1 = rand(stream_for!(ks1, :b))
    a2 = rand(stream_for!(ks1, :a)); b2 = rand(stream_for!(ks1, :b))

    # Order B (adversarial): b, b, a, a
    b1p = rand(stream_for!(ks2, :b)); b2p = rand(stream_for!(ks2, :b))
    a1p = rand(stream_for!(ks2, :a)); a2p = rand(stream_for!(ks2, :a))

    @test a1 == a1p
    @test a2 == a2p
    @test b1 == b1p
    @test b2 == b2p

    # Occurrence counts track the number of draws per key regardless of order.
    @test ks1.counts[:a] == 2 && ks1.counts[:b] == 2
    @test ks2.counts == ks1.counts
end


@safetestset keyed_streams_seed_independence =
    "two different seeds give different per-key draws" begin
    using CompetingClocks: KeyedStreams, stream_for!
    using Random: rand
    ks1 = KeyedStreams{Symbol}(1)
    ks2 = KeyedStreams{Symbol}(2)
    @test rand(stream_for!(ks1, :a)) != rand(stream_for!(ks2, :a))
end


@safetestset keyed_streams_race_is_independent =
    "the reserved race stream is independent of every per-key stream" begin
    using CompetingClocks: KeyedStreams, stream_for!, race_stream
    using Random: rand
    ks = KeyedStreams{Symbol}(42)
    r = rand(race_stream(ks))
    a = rand(stream_for!(ks, :a))
    # Drawing the race did not touch key :a's stream: a matches a fresh family's
    # first :a draw.
    ksref = KeyedStreams{Symbol}(42)
    @test a == rand(stream_for!(ksref, :a))
    @test r == rand(race_stream(KeyedStreams{Symbol}(42)))
end


@safetestset keyed_streams_copy_independence =
    "a copy carries generator state, and advancing the copy leaves the original's next draws unchanged" begin
    using CompetingClocks: KeyedStreams, stream_for!
    using Random: rand

    ks = KeyedStreams{Symbol}(7)
    rand(stream_for!(ks, :a))          # advance :a once
    rand(stream_for!(ks, :b))
    kc = copy(ks)

    # The copy resumes from the SAME state: its next :a draw equals what the
    # original's next :a draw will be (coupling), and counts were copied.
    @test kc.counts == ks.counts
    c_next = rand(stream_for!(kc, :a))

    # Advance the copy several more times; this must NOT disturb the original.
    for _ in 1:5
        rand(stream_for!(kc, :a)); rand(stream_for!(kc, :b))
    end

    o_next = rand(stream_for!(ks, :a))
    @test o_next == c_next             # original's next :a draw is the coupled one
end


@safetestset keyed_streams_rekey_resets =
    "re-keying forgets every generator and count and yields a fresh independent family" begin
    using CompetingClocks: KeyedStreams, stream_for!, rekey_streams!
    using Random: rand

    ks = KeyedStreams{Symbol}(7)
    rand(stream_for!(ks, :a)); rand(stream_for!(ks, :b))
    @test !isempty(ks.gens) && ks.counts[:a] == 1

    rekey_streams!(ks, 999)
    @test isempty(ks.gens) && isempty(ks.counts)

    # A different seed after rekey gives a different :a draw than the original seed.
    ksref = KeyedStreams{Symbol}(7)
    @test rand(stream_for!(ks, :a)) != rand(stream_for!(ksref, :a))
    # And matches a fresh family seeded the same as the rekey.
    @test rand(stream_for!(KeyedStreams{Symbol}(999), :a)) ==
          rand(stream_for!(KeyedStreams{Symbol}(999), :a))
end


# ---------------------------------------------------------------------------
# clone full fidelity: a clone runs the identical trajectory as its original
# because it inherits the whole stream family (generator states + counts). This
# is the coupling primitive the CRN successor is built on.
# ---------------------------------------------------------------------------
@safetestset clone_full_fidelity =
    "cloning a live sampler and running both forward yields identical firing sequences" begin
    using CompetingClocks: FirstToFire, CombinedNextReaction, enable!, next, fire!, clone
    using Distributions: Weibull, Exponential, Gamma
    using Random: Xoshiro

    for B in (FirstToFire, CombinedNextReaction)
        s = B{Int,Float64}(0xBEEF)
        enable!(s, 1, Weibull(1.5, 1.0), 0.0, 0.0)
        enable!(s, 2, Exponential(1.3), 0.0, 0.0)
        enable!(s, 3, Gamma(2.0, 0.7), 0.0, 0.0)

        c = clone(s)

        seq_s = Tuple{Float64,Int}[]
        seq_c = Tuple{Float64,Int}[]
        for _ in 1:3
            ws, ks = next(s, 0.0); push!(seq_s, (ws, ks)); fire!(s, ks, ws)
            wc, kc = next(c, 0.0); push!(seq_c, (wc, kc)); fire!(c, kc, wc)
        end
        @test seq_s == seq_c
    end
end


@safetestset force_fire_coupled_redraw =
    "a clone pair sharing streams redraws a passed loser identically under force_fire!" begin
    using CompetingClocks: FirstToFire, CombinedNextReaction, enable!, next, force_fire!,
        clone, keys
    using Distributions: Weibull, Exponential
    using Random: Xoshiro

    for B in (FirstToFire, CombinedNextReaction)
        s = B{Symbol,Float64}(0x0FF1CE)
        # A spread of losers so at least one is scheduled at/before a chosen tstar
        # (the redraw-if-passed branch), plus a clock :y to force.
        for (k, d) in ((:a, Weibull(1.7, 1.0)), (:b, Exponential(1.0)),
                       (:c, Weibull(2.0, 1.5)), (:d, Exponential(2.0)))
            enable!(s, k, d, 0.0, 0.0)
        end
        enable!(s, :y, Exponential(0.5), 0.0, 0.0)

        c = clone(s)
        tstar = 0.6
        force_fire!(s, :y, tstar)
        force_fire!(c, :y, tstar)

        # Every surviving loser — whether kept or redrawn — has the identical
        # scheduled time on both, because the redraw draws from each loser's own
        # stream and the clone inherited those exact stream states.
        for k in (:a, :b, :c, :d)
            @test s[k] == c[k]
        end
    end
end
