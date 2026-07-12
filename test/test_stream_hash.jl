using SafeTestsets

# Unit tests for the `stream_hash` seam (phase OB-3a). The design finding these
# certify: every per-key stream-seed derivation routes through ONE overloadable
# function, so (a) the default behavior is bit-for-bit what `hash((seed, key))`
# always produced, and (b) a caller can give an alternate key REPRESENTATION the
# same stream identity as its content tuple by overloading the seam for its own
# key type — without touching any KeyedStreams internals.


@safetestset stream_hash_default_matches_tuple_hash =
    "the default stream_hash is bit-for-bit hash((seed, key)), so pre-seam stream derivations are unchanged" begin
    using CompetingClocks: KeyedStreams, stream_for!, stream_hash
    using Random: Xoshiro, rand

    seed = UInt64(0x9E3779B97F4A7C15)
    for key in [(:a, 3), (:b,), (:machine, 7, 2.5)]
        @test stream_hash(seed, key) === hash((seed, key))
    end

    # A live family's first draw for a key equals a generator hand-seeded from
    # the historical derivation: the seam changed no bytes on the default path.
    ks = KeyedStreams{Tuple}(seed)
    manual = Xoshiro(hash((seed, (:a, 3))))
    @test rand(stream_for!(ks, (:a, 3))) == rand(manual)
end


@safetestset stream_hash_overload_changes_derivation =
    "overloading stream_hash for a key type re-derives that type's streams while leaving other keys' streams alone" begin
    using CompetingClocks: KeyedStreams, stream_for!
    import CompetingClocks
    using Random: Xoshiro, rand

    # A struct key standing in for its content tuple. The struct's default hash
    # is NOT its tuple's hash, so without the overload the two representations
    # would draw different streams from the same family seed. The definition is
    # QUALIFIED (CompetingClocks.stream_hash): the @testset body is a soft
    # local scope, where an unqualified `stream_hash(...) = ...` would silently
    # create a local function instead of extending the seam.
    struct StructKey
        name::Symbol
        idx::Int
    end
    CompetingClocks.stream_hash(seed::UInt64, key::StructKey) =
        hash((seed, (key.name, key.idx)))

    seed = UInt64(42)
    @test hash(StructKey(:a, 3)) != hash((:a, 3))   # the motivating mismatch

    # KeyedStreams keyed by the struct produce the SAME first draw as a family
    # keyed by the content tuple: the overload is the whole cross-representation
    # stream-identity mechanism. invokelatest: the overload was defined inside
    # this testset's single top-level evaluation, so calls must cross the world
    # age it created (a package overloading at its own top level needs nothing).
    ks_struct = KeyedStreams{StructKey}(seed)
    ks_tuple = KeyedStreams{Tuple}(seed)
    draw_struct = Base.invokelatest(() -> rand(stream_for!(ks_struct, StructKey(:a, 3))))
    draw_tuple = Base.invokelatest(() -> rand(stream_for!(ks_tuple, (:a, 3))))
    @test draw_struct == draw_tuple
end


@safetestset stream_hash_seam_reaches_pinned_keys =
    "a pinned key's replay-from-seed derivation goes through the same stream_hash seam as first use" begin
    using CompetingClocks: KeyedStreams, stream_for!, pin_stream!, rekey_streams!
    using Random: rand

    # Pin a key, draw, rekey the family: the pinned key must REPLAY its original
    # stream, which only happens if the pinned re-derivation uses the same
    # (stream_hash-routed) derivation as the first one.
    ks = KeyedStreams{Tuple}(UInt64(7))
    pin_stream!(ks, (:init,))
    first_draw = rand(stream_for!(ks, (:init,)))
    rekey_streams!(ks, UInt64(999))
    @test rand(stream_for!(ks, (:init,))) == first_draw
end
