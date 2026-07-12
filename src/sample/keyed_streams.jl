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
holds each key's occurrence count (how many times it has drawn). `pinned` maps
a key to the family seed its stream derives from FOREVER (see
[`pin_stream!`](@ref)): a pinned key's draws are immune to
[`rekey_streams!`](@ref), which is how a caller keeps one reserved stream (e.g.
an initialization stream) shared across otherwise-decoupled clones.

!!! warning "Keys must hash by content"
    Stream seeding uses [`stream_hash`](@ref)`(seed, key)` (by default
    `hash((seed, key))`), so seeded reproducibility requires every key
    component to have a CONTENT-based `Base.hash`. Symbols, numbers,
    strings, and tuples of them qualify. An `@enum` does NOT: it falls back to
    `objectid`, which differs across processes (and even across compile
    options), silently breaking same-seed reproducibility and cross-module key
    identity. If a clock key contains an enum `E`, define
    `Base.hash(x::E, h::UInt) = hash(Symbol(x), h)` next to the enum. A caller
    whose key type needs a DIFFERENT content identity (e.g. an event struct
    standing in for its field tuple) overloads `stream_hash` for that type;
    the obligation that the hashed content be content-hashable is the same.
"""
mutable struct KeyedStreams{K}
    seed::UInt64
    gens::Dict{K,Xoshiro}
    counts::Dict{K,Int}
    race_gen::Xoshiro
    pinned::Dict{K,UInt64}
end

function KeyedStreams{K}(seed=_DEFAULT_STREAM_SEED) where {K}
    s = UInt64(seed)
    return KeyedStreams{K}(
        s, Dict{K,Xoshiro}(), Dict{K,Int}(), Xoshiro(hash((s, _RACE_SENTINEL))),
        Dict{K,UInt64}(),
    )
end

# The seed a key's generator derives from: the pinned seed if the key is
# pinned, else the current family seed.
_stream_seed(ks::KeyedStreams{K}, key::K) where {K} = get(ks.pinned, key, ks.seed)

"""
    stream_hash(seed::UInt64, key) -> UInt64

The canonical per-key stream-seed derivation: the number that seeds `key`'s
`Xoshiro` generator in a [`KeyedStreams`](@ref) family with family seed `seed`.
The default is `hash((seed, key))`, and every (seed, key) hash site in this
package routes through this ONE function, so it is the overloadable seam that
defines a key's stream identity.

Why a seam: a framework may address the same clock by more than one key
REPRESENTATION — for example, an event struct instance instead of the tuple of
its type name and fields. A struct's default `hash` is not its content tuple's
`hash`, so without an overload the two representations would draw DIFFERENT
streams from the same seed. The framework restores cross-representation stream
identity by defining, for its own key type,

```julia
CompetingClocks.stream_hash(seed::UInt64, key::MyKey) = hash((seed, content_tuple(key)))
```

so both representations of one clock derive the same generator.

The content-hash obligation stays with the caller either way: whatever value
the overload hashes must itself hash by CONTENT (see the warning on
[`KeyedStreams`](@ref) — an `@enum` field hashes by `objectid` and silently
breaks same-seed reproducibility unless given a content `Base.hash`). The
reserved race stream does not route through this seam; its sentinel key is
package-internal and never has an alternate representation.
"""
stream_hash(seed::UInt64, key) = hash((seed, key))

"""
    stream_for!(ks, key) -> AbstractRNG

Return `key`'s generator, creating it (seeded from `(seed, key)`) on first use,
and bump `key`'s occurrence count. The caller draws from the returned generator
directly — `rand(gen, dist)` may consume several native words, but they all come
from THIS key's stream, so the key's future is independent of every other key's.
A pinned key (see [`pin_stream!`](@ref)) derives from its pinned seed instead of
the current family seed.
"""
function stream_for!(ks::KeyedStreams{K}, key::K) where {K}
    gen = get!(ks.gens, key) do
        # The (seed, key) hash routes through the `stream_hash` seam so a caller
        # can pin a key's stream identity to its CONTENT rather than its
        # representation (see the stream_hash docstring).
        Xoshiro(stream_hash(_stream_seed(ks, key), key))
    end
    ks.counts[key] = get(ks.counts, key, 0) + 1
    return gen
end

"""
    pin_stream!(ks, key) -> ks

Pin `key`'s seed derivation to the CURRENT family seed. From now on, `key`'s
generator always derives from `hash((pinned_seed, key))`, even after any number
of [`rekey_streams!`](@ref) calls: rekeying still discards `key`'s live
generator state (like every other key's), but the key re-derives from its
pinned seed, so it REPLAYS its original stream from the start instead of moving
to the new family. `copy` carries pins, so a copied family's pinned keys keep
the same derivation.

This is the carve-out a simulation framework uses for a reserved
initialization stream: decoupling a clone with a new seed must give the clone
independent futures WITHOUT re-deriving the initial condition's draws —
branching happens after time zero, so clones share the time-zero randomness.
Pinning the same key again re-records the pin at the family seed current at
that moment.

`pin_stream!` is generic over keys; the caller decides which key is reserved.
"""
function pin_stream!(ks::KeyedStreams{K}, key::K) where {K}
    ks.pinned[key] = ks.seed
    return ks
end

"""
    is_pinned(ks, key) -> Bool

Whether `key`'s seed derivation is pinned (see [`pin_stream!`](@ref)).
"""
is_pinned(ks::KeyedStreams{K}, key::K) where {K} = haskey(ks.pinned, key)

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

A PINNED key (see [`pin_stream!`](@ref)) keeps its seed derivation: its live
generator state and count are discarded like every other key's, but its next
draw re-derives from the pinned seed, replaying the pinned key's ORIGINAL
stream from the start rather than moving to the new family. Per-key streams
replay from a seed alone, so preserving only the derivation is enough to keep
the pinned key's randomness shared between a clone and its original.
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
against — so both are copied here in one place. Pins are copied too, so the
copy's pinned keys keep the same seed derivation across its own rekeys.
"""
function Base.copy(ks::KeyedStreams{K}) where {K}
    gens = Dict{K,Xoshiro}(k => copy(v) for (k, v) in ks.gens)
    return KeyedStreams{K}(ks.seed, gens, copy(ks.counts), copy(ks.race_gen), copy(ks.pinned))
end
