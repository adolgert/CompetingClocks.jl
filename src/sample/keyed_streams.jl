# ---------------------------------------------------------------------------
# Keyed random streams: the sampler owns randomness, addressed by clock.
#
# The randomness-ownership inversion (milestone 4). No sampler verb takes an
# RNG anymore. Instead every sampler holds one `KeyedStreams`, and a draw
# belongs to a CLOCK (and an occurrence), not to a position in a global call
# sequence. Two runs built from the SAME seed therefore consume identical
# randomness per clock even when the event ORDER differs, because each key's
# generator is seeded once from the key and advances only when that key draws.
# This per-clock addressing is what makes clone coupling and common random
# numbers compositional: coupling a θ and a θ+h run, or replaying one event
# order against another, needs nothing but the same seed.
#
# Per-key seeding uses `Xoshiro(hash((seed, key)))`: the hash runs ONCE per key
# to pick that key's seed, never once per uniform. A real generator (Xoshiro)
# then produces the key's stream. This is standard practice; hashing per uniform
# to synthesize randomness is not.
#
# Reserved stream for race-owned draws: some samplers (DirectCall, MultipleDirect,
# RSSA, PSSACR, Petri) draw for the RACE as a whole — the time to the next event
# and which clock wins — rather than for one clock. Those draws have no natural
# clock key, so instead of synthesizing a reserved value of the (arbitrary) key
# type K, `KeyedStreams` carries a dedicated `race_gen` field beside the keyed
# table. `race_stream(ks)` hands it back. Keeping it inside `KeyedStreams` means
# `copy` and `rekey_streams!` handle the race generator in the same place as the
# per-key ones, so no caller can forget to copy or reseed it — that omission is
# exactly the silent-coupling bug class this design exists to prevent.
# ---------------------------------------------------------------------------

using Random: Xoshiro, AbstractRNG

# Default seed used by a sampler constructed without an explicit one. A context
# immediately rekeys the sampler it builds from its own rng, and tests that care
# about a particular stream pass a seed explicitly; the default only keeps a bare
# `Sampler{K,T}()` deterministic.
const _DEFAULT_STREAM_SEED = 0x243F6A8885A308D3

# Distinguishes the race generator's seed from any per-key seed. A key of type K
# can never collide with this Symbol inside `hash`, so the race stream is
# independent of every per-key stream.
const _RACE_SENTINEL = :__competingclocks_race_stream__

"""
    KeyedStreams{K}(seed)

Per-clock random streams owned by a sampler. `seed` selects the whole family of
streams: key `k`'s generator is `Xoshiro(hash((seed, k)))` and the race
generator is `Xoshiro(hash((seed, :__competingclocks_race_stream__)))`. `gens`
holds the live per-key generators (created lazily on first draw) and `counts`
holds each key's occurrence count (how many times it has drawn).

!!! warning "Keys must hash by content"
    Stream seeding uses `hash((seed, key))`, so seeded reproducibility requires
    every key component to have a CONTENT-based `Base.hash`. Symbols, numbers,
    strings, and tuples of them qualify. An `@enum` does NOT: it falls back to
    `objectid`, which differs across processes (and even across compile
    options), silently breaking same-seed reproducibility and cross-module key
    identity. If a clock key contains an enum `E`, define
    `Base.hash(x::E, h::UInt) = hash(Symbol(x), h)` next to the enum.
"""
mutable struct KeyedStreams{K}
    seed::UInt64
    gens::Dict{K,Xoshiro}
    counts::Dict{K,Int}
    race_gen::Xoshiro
end

function KeyedStreams{K}(seed=_DEFAULT_STREAM_SEED) where {K}
    s = UInt64(seed)
    return KeyedStreams{K}(s, Dict{K,Xoshiro}(), Dict{K,Int}(), Xoshiro(hash((s, _RACE_SENTINEL))))
end

"""
    stream_for!(ks, key) -> AbstractRNG

Return `key`'s generator, creating it (seeded from `(seed, key)`) on first use,
and bump `key`'s occurrence count. The caller draws from the returned generator
directly — `rand(gen, dist)` may consume several native words, but they all come
from THIS key's stream, so the key's future is independent of every other key's.
"""
function stream_for!(ks::KeyedStreams{K}, key::K) where {K}
    gen = get!(ks.gens, key) do
        Xoshiro(hash((ks.seed, key)))
    end
    ks.counts[key] = get(ks.counts, key, 0) + 1
    return gen
end

"""
    race_stream(ks) -> AbstractRNG

The reserved generator for draws that belong to the race rather than a clock
(next-event time, winner selection). Independent of every per-key stream.
"""
race_stream(ks::KeyedStreams) = ks.race_gen

"""
    rekey_streams!(ks, seed)

Re-seed the whole family to `seed` and forget every live per-key generator and
count. After this, the streams produce a fresh, independent trajectory of draws;
this is the mechanism that decouples a clone from its original and that gives
each split copy its own randomness.
"""
function rekey_streams!(ks::KeyedStreams{K}, seed) where {K}
    ks.seed = UInt64(seed)
    empty!(ks.gens)
    empty!(ks.counts)
    ks.race_gen = Xoshiro(hash((ks.seed, _RACE_SENTINEL)))
    return ks
end

"""
    copy(ks)

A deep copy carrying the GENERATOR STATES and counts, not just the seed.
Advancing the copy leaves the original's next draws untouched and vice versa.
Forgetting to copy a generator's state (or a count) would silently couple the
copy's future to the original's — the known bug class this whole design guards
against — so both are copied here in one place.
"""
function Base.copy(ks::KeyedStreams{K}) where {K}
    gens = Dict{K,Xoshiro}(k => copy(v) for (k, v) in ks.gens)
    return KeyedStreams{K}(ks.seed, gens, copy(ks.counts), copy(ks.race_gen))
end
