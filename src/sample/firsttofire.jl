using DataStructures: MutableBinaryMinHeap, extract_all!, update!

# Helper struct for FirstToFire to be able to resample.
struct FTFEntry{T,D}
    handle::Int
    distribution::D
    te::T
end


"""
    FirstToFire{KeyType,TimeType}(seed; coupling=:carry)

This sampler is often the fastest for non-exponential distributions.
When a clock is first enabled, this sampler asks the clock when it would
fire and saves that time in a sorted heap of future times. Then it works
through the heap, one by one. When a clock is disabled, its future firing time
is removed from the list. There is no memory of previous firing times.

This uses a `DataStructures.MutableBinaryMinHeap` which is a Fibonacci
heap. It has been tested against many other heaps and rarely loses by more
than a few percent, so we are sticking with it.

The `coupling` keyword (default `:carry`) fixes, at construction, how
[`reenable!`](@ref) realizes a mid-flight distribution change: `:carry` rebases
the stored firing time through the change deterministically (consuming no
randomness — the coupling pathwise/IPA derivatives need), while `:redraw` draws
the remaining lifetime fresh, conditioned on age. Both agree in law. Read it
back with [`coupling`](@ref).
"""
mutable struct FirstToFire{K,T} <: SSA{K,T}
    firing_queue::MutableBinaryMinHeap{OrderedSample{K,T}}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{K,FTFEntry{T,UnivariateDistribution}}
    streams::KeyedStreams{K}
    coupling::Symbol
end


function FirstToFire{K,T}(seed=_DEFAULT_STREAM_SEED; coupling::Symbol=:carry) where {K,T}
    validate_coupling(FirstToFire, coupling)
    heap = MutableBinaryMinHeap{OrderedSample{K,T}}()
    state = Dict{K,FTFEntry{T,UnivariateDistribution}}()
    FirstToFire{K,T}(heap, state, KeyedStreams{K}(seed), coupling)
end

coupling(ftf::FirstToFire) = ftf.coupling

# A full-state clone copies the firing queue, the entry table, AND the keyed
# streams (generator states, counts). Running the clone and the original forward
# yields identical firing sequences — every clock was drawn once at enable from
# its own stream and the clone inherits those exact stream states.
function clone(ftf::FirstToFire{K,T}) where {K,T}
    c = FirstToFire{K,T}(ftf.streams.seed; coupling=ftf.coupling)
    c.firing_queue = deepcopy(ftf.firing_queue)
    copy!(c.transition_entry, ftf.transition_entry)
    c.streams = copy(ftf.streams)
    return c
end

similar_sampler(ftf::FirstToFire{K,T}) where {K,T} =
    FirstToFire{K,T}(ftf.streams.seed; coupling=ftf.coupling)

rekey_streams!(ftf::FirstToFire, seed) = (rekey_streams!(ftf.streams, seed); ftf)


function reset!(propagator::FirstToFire{K,T}) where {K,T}
    # Drain by pop! rather than empty!: DataStructures gained
    # empty!(::MutableBinaryMinHeap) only in 0.19, and the compat range admits
    # 0.18 so this package can co-resolve with Gen.jl's 0.18 cap.
    while !isempty(propagator.firing_queue)
        pop!(propagator.firing_queue)
    end
    empty!(propagator.transition_entry)
end

function copy_clocks!(dst::FirstToFire{K,T}, src::FirstToFire{K,T}) where {K,T}
    dst.firing_queue = deepcopy(src.firing_queue)
    copy!(dst.transition_entry, src.transition_entry)
    dst.streams = copy(src.streams)
    # copy_clocks! promises a full replacement of the destination's state, so
    # the destination must re-evaluate clocks the same way the source would.
    dst.coupling = src.coupling
    return dst
end


function jitter!(propagator::FirstToFire{K,T}, when::T) where {K,T}
    for (clock, entry) in propagator.transition_entry
        te = entry.te
        distribution = entry.distribution
        gen = stream_for!(propagator.streams, clock)
        if te < when
            when_fire = te + rand(gen, truncated(distribution, when - te, typemax(T)))
        else
            when_fire = te + rand(gen, distribution)
        end
        update!(propagator.firing_queue, entry.handle, OrderedSample{K,T}(clock, when_fire))
    end
end


# Finds the next one without removing it from the queue.
function next(propagator::FirstToFire{K,T}, when::T) where {K,T}
    least = if !isempty(propagator.firing_queue)
        first(propagator.firing_queue)
    else
        OrderedSample(nothing, typemax(T))
    end
    (least.time, least.key)
end


function enable!(
    propagator::FirstToFire{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T) where {K,T}

    gen = stream_for!(propagator.streams, clock)
    if te < when
        when_fire = te + rand(gen, truncated(distribution, when - te, typemax(T)))
    else
        when_fire = te + rand(gen, distribution)
    end
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock].handle
        update!(propagator.firing_queue, heap_handle, OrderedSample{K,T}(clock, when_fire))
    else
        heap_handle = push!(propagator.firing_queue, OrderedSample{K,T}(clock, when_fire))
    end
    propagator.transition_entry[clock] = FTFEntry{T,UnivariateDistribution}(
        heap_handle,
        distribution,
        te
    )
end


# fire! is the SSA interface fallback (forwards to disable!): FirstToFire
# holds a clock's drawn firing time only while the clock is enabled and
# discards it at disable!, so it retains no residual draw randomness and
# firing and disabling act identically on its state.


function disable!(propagator::FirstToFire{K,T}, clock::K, when::T) where {K,T}
    entry = get(propagator.transition_entry, clock, nothing)
    entry === nothing && throw(KeyError(clock))
    delete!(propagator.firing_queue, entry.handle)
    delete!(propagator.transition_entry, clock)
end

"""
    getindex(sampler::FirstToFire{K,T}, clock::K)

For the `FirstToFire` sampler, returns the stored firing time associated to the clock.
"""
function Base.getindex(propagator::FirstToFire{K,T}, clock::K) where {K,T}
    if haskey(propagator.transition_entry, clock)
        heap_handle = propagator.transition_entry[clock].handle
        return getfield(propagator.firing_queue[heap_handle], :time)
    else
        throw(KeyError(clock))
    end
end

function Base.keys(propagator::FirstToFire)
    return collect(keys(propagator.transition_entry))
end

function Base.length(propagator::FirstToFire)
    return length(propagator.transition_entry)
end


Base.haskey(propagator::FirstToFire{K,T}, clock::K) where {K,T} = haskey(propagator.transition_entry, clock)
Base.haskey(propagator::FirstToFire{K,T}, clock) where {K,T} = false

enabled(propagator::FirstToFire) = keys(propagator.transition_entry)
isenabled(propagator::FirstToFire{K,T}, clock::K) where {K,T} = haskey(propagator.transition_entry, clock)


# --- estimator-facing verbs -------------------------------------------------
# FirstToFire is a scheduling backend: it draws once at enable and stores the
# absolute firing time in a heap. Every entry in `transition_entry` is enabled
# (disable! deletes it), so its enabling times feed `enabled_ages` directly, and
# it supports forcing with the keep-if-later / redraw-if-passed repair. It does
# NOT retain the draw's survival uniform — the heap holds only the firing time —
# so `retained_draw` stays unsupported.

supports_enabled_ages(::Type{<:FirstToFire}) = true
supports_force(::Type{<:FirstToFire}) = true
# FirstToFire stores only the absolute firing time, not a survival uniform, so it
# cannot answer retained_draw — but it CAN carry: the retained log-uniform of the
# total lifetime is reconstructible from the stored time as logccdf(dist, sched -
# te), which is all the rebasing map needs. So it supports the deterministic-carry
# coupling even though supports_retained_draw stays false.
supports_carry(::Type{<:FirstToFire}) = true

enabling_times(propagator::FirstToFire) =
    ((clock, entry.te) for (clock, entry) in propagator.transition_entry)


"""
    reenable!(propagator::FirstToFire, clock, distribution, te, when)

Re-evaluate `clock`'s distribution mid-flight. Which pathwise coupling realizes
the change is the sampler's construction-time [`coupling`](@ref) field.

`:redraw` is the plain [`enable!`](@ref)-on-enabled-key path: draw a fresh firing
time from the clock's own keyed stream, conditioned on the current age.

`:carry` reconstructs the retained draw from the stored firing time and rebases
it — the five-liner of the design contract, consuming no randomness. With `a` the
age at the change and `sched_old` the stored firing time, the retained
log-uniform of the TOTAL lifetime is `logu = logccdf(dist_old, sched_old − te)`
(FirstToFire keeps no uniform, but the stored time determines it); the new firing
time then solves conditional-survival matching:

    lsurv_new = logu − logccdf(dist_old, a) + logccdf(dist_new, a_new)
    sched_new = te + invlogccdf(dist_new, lsurv_new)

where `a` is measured from the OLD enabling time and `a_new` from the (possibly
re-anchored) new `te`. For exponentials this reduces exactly to Gibson–Bruck
rate-ratio scaling of the remaining time.
"""
function reenable!(
    propagator::FirstToFire{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T) where {K,T}
    haskey(propagator, clock) || throw(ArgumentError(
        "reenable! needs clock $clock currently enabled; use enable! to start a " *
        "clock that is not enabled."))
    # The constructor validated the field, so only the two couplings reach here.
    if propagator.coupling === :carry
        _reenable_carry!(propagator, clock, distribution, te, when)
    else
        enable!(propagator, clock, distribution, te, when)
    end
    nothing
end


function _reenable_carry!(
    propagator::FirstToFire{K,T}, clock::K, distribution::UnivariateDistribution,
    te::T, when::T) where {K,T}
    entry = propagator.transition_entry[clock]
    dist_old = entry.distribution
    te_old = entry.te
    sched_old = getfield(propagator.firing_queue[entry.handle], :time)
    a_old = when - te_old               # age at the change under the old anchor
    a_new = when - te                   # age at the change under the new anchor
    # The retained log-uniform of the total lifetime, recovered from the stored
    # firing time: sched_old == te_old + invlogccdf(dist_old, logu).
    logu = logccdf(dist_old, sched_old - te_old)
    lsurv_new = logu - logccdf(dist_old, a_old) + logccdf(distribution, a_new)
    sched_new = te + invlogccdf(distribution, lsurv_new)
    update!(propagator.firing_queue, entry.handle, OrderedSample{K,T}(clock, sched_new))
    propagator.transition_entry[clock] = FTFEntry{T,UnivariateDistribution}(
        entry.handle, distribution, te
    )
    nothing
end

"""
    force_fire!(propagator::FirstToFire, key, tstar)

Fire `key` at `tstar`. The fired clock is removed (mirroring `fire!`, which
forwards to `disable!`). Each surviving clock is repaired: keep-if-later leaves a
scheduled time already past `tstar` untouched — conditioned on the race being
decided at `tstar` an unconditional draw that exceeds `tstar` is a sample from
the survival-conditioned law — while redraw-if-passed replaces a scheduled time
`≤ tstar` (a promise the force overran) with a fresh draw from the clock's
lifetime truncated past `tstar`, drawn from that loser's OWN keyed stream. Using
the loser's stream is what makes the redraw identical across a coupled clone pair.
"""
function force_fire!(propagator::FirstToFire{K,T}, key::K, tstar::T) where {K,T}
    entry = get(propagator.transition_entry, key, nothing)
    entry === nothing && throw(KeyError(key))
    delete!(propagator.firing_queue, entry.handle)
    delete!(propagator.transition_entry, key)
    for (k, e) in propagator.transition_entry
        scheduled = getfield(propagator.firing_queue[e.handle], :time)
        scheduled > tstar && continue  # keep-if-later
        # A passed promise: redraw from the lifetime truncated past tstar. The
        # age tstar - te is nonnegative here because scheduled ≥ te and
        # scheduled ≤ tstar.
        shifted = truncated(e.distribution, tstar - e.te, typemax(T))
        newtime = e.te + rand(stream_for!(propagator.streams, k), shifted)
        update!(propagator.firing_queue, e.handle, OrderedSample{K,T}(k, newtime))
    end
    nothing
end
