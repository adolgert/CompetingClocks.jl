# ---------------------------------------------------------------------------
# Estimator-facing sampler verbs and the capability traits that gate them.
#
# A derivative estimator (score, pathwise/IPA, weak-derivative branching)
# consumes more than the four simulation verbs enable!/disable!/next/fire!. It
# needs to ask a live sampler which clocks are enabled and how aged they are
# (`enabled_ages`), to read the survival-space uniform behind a scheduled draw
# (`retained_draw`), and to fire a CHOSEN clock at a CHOSEN time regardless of
# the race (`force_fire!`). Not every sampler can answer every question: the
# exponential-only samplers keep no per-clock enabling-time table, so age is
# meaningless for them; only the scheduling backend that retains draw
# randomness can hand back a retained uniform. These traits let an estimator
# check what a sampler supports at construction, and let a verb throw a clear
# ArgumentError naming the missing trait rather than a MethodError deep inside a
# backend.
#
# The traits follow the Holy-traits style already used by `sampling_space`
# (src/sample/nrtransition.jl): a method on the concrete type decides the
# capability, an instance method forwards to it, and the default is `false` so a
# newly added sampler opts out until it deliberately opts in.
# ---------------------------------------------------------------------------

"""
    supports_enabled_ages(sampler) -> Bool

Whether the sampler can answer [`enabled_ages`](@ref): whether it stores a
per-clock enabling-time table. Age (`when - te`) is only meaningful for a
sampler that carries lifetime distributions with memory; the exponential-only
samplers (`DirectCall`, `MultipleDirect`, `RSSA`, `PSSACR`) and `Petri` re-read
a memoryless rate from the state each step and keep no enabling times, so they
leave this `false`.
"""
supports_enabled_ages(::Type{<:SSA}) = false
supports_enabled_ages(s::SSA) = supports_enabled_ages(typeof(s))

"""
    supports_force(sampler) -> Bool

Whether the sampler implements [`force_fire!`](@ref), firing a chosen clock at a
chosen time and re-conditioning the losers on survival past that time.
"""
supports_force(::Type{<:SSA}) = false
supports_force(s::SSA) = supports_force(typeof(s))

"""
    supports_retained_draw(sampler) -> Bool

Whether the sampler implements [`retained_draw`](@ref), returning the
survival-space uniform behind a clock's current tentative firing time. Only a
scheduling backend that RETAINS draw randomness across state changes can, which
in this package is `CombinedNextReaction` alone: `FirstToFire` discards its draw
at every disable, and `FirstReaction` redraws from scratch on every `next`.
"""
supports_retained_draw(::Type{<:SSA}) = false
supports_retained_draw(s::SSA) = supports_retained_draw(typeof(s))

"""
    supports_carry(sampler) -> Bool

Whether the sampler can implement the DETERMINISTIC-CARRY re-evaluation coupling
of [`reenable!`](@ref): when a still-enabled clock's distribution changes
mid-flight, map its retained draw through the change by matching conditional
survival, consuming NO fresh randomness. This is the only coupling under which a
firing time moves CONTINUOUSLY in a distribution parameter, so it is the one
pathwise (infinitesimal perturbation analysis) derivatives need.

A sampler can carry only if it retains enough of the in-flight draw to rebase it:
`CombinedNextReaction` keeps the survival uniform outright, and `FirstToFire`
reconstructs it from the stored absolute firing time (its retained-draw identity
`logu = logccdf(dist_old, sched - te)`). The exponential-only samplers keep no
in-flight draw to carry (they re-read a memoryless rate each step), and
`FirstReaction` redraws every clock at every `next` so it has no retained
schedule to move — both leave this `false`, and constructing them with
`coupling=:carry` throws an `ArgumentError` naming this trait (see
[`validate_coupling`](@ref)).
"""
supports_carry(::Type{<:SSA}) = false
supports_carry(s::SSA) = supports_carry(typeof(s))


"""
    validate_coupling(SamplerType, coupling::Symbol) -> Symbol

Check a re-evaluation `coupling` requested at sampler construction: it must be
`:carry` or `:redraw`, and `:carry` additionally requires
[`supports_carry`](@ref)`(SamplerType) == true`. Throws an `ArgumentError`
otherwise; returns `coupling` so a constructor can validate and store in one
expression. Every sampler constructor (and every `SamplerSpec` that forwards a
coupling) routes through this ONE helper so the validation rules live in a
single place and a bad request fails at construction, not at the first
[`reenable!`](@ref).
"""
function validate_coupling(::Type{S}, coupling::Symbol) where {S<:SSA}
    coupling in (:carry, :redraw) || throw(ArgumentError(
        "coupling must be :carry or :redraw, got :$coupling."))
    if coupling === :carry && !supports_carry(S)
        throw(ArgumentError(
            "coupling=:carry is not supported by $S; supports_carry(::Type{$S}) " *
            "is false. This sampler keeps no in-flight draw to carry through a " *
            "distribution change; construct it with coupling=:redraw, or choose " *
            "a sampler with supports_carry == true."))
    end
    return coupling
end


"""
    coupling(sampler) -> Symbol

The re-evaluation coupling this sampler applies when [`reenable!`](@ref)
changes a still-enabled clock's distribution: `:carry` (map the retained draw
through the change deterministically) or `:redraw` (draw the remaining lifetime
fresh, conditioned on age). The coupling is a CONSTRUCTION-TIME property of the
sampler, chosen by the `coupling` keyword of the carry-capable samplers
(`CombinedNextReaction`, `FirstToFire`, default `:carry`). Every other sampler
retains no in-flight draw, so redraw is its only possible behavior and this
generic method answers `:redraw` for it. A `SamplingContext` forwards to its
inner sampler. Downstream recorders read a run's coupling from this accessor.
"""
coupling(::SSA) = :redraw


"""
    enabling_times(sampler)

Iterate `(key, te)` pairs for every clock the sampler currently holds enabled,
where `te` is the clock's enabling time (the zero of its lifetime distribution,
in absolute time). This is the one per-sampler accessor [`enabled_ages`](@ref)
is written against; a sampler supplies it wherever the information already lives
in its clock table. Only defined for samplers with
`supports_enabled_ages == true`.
"""
function enabling_times(s::SSA)
    throw(ArgumentError(
        "enabling_times is not available for $(typeof(s)); " *
        "supports_enabled_ages(::$(typeof(s))) is $(supports_enabled_ages(s)). " *
        "Only samplers that keep a per-clock enabling-time table implement it."))
end


"""
    enabled_ages(sampler, when) -> Vector{Tuple{K,Float64}}

Every currently-enabled clock paired with its age `when - te`, SORTED BY KEY so
the order is deterministic across clones and backends. The deterministic order
is load-bearing: a forced-selection probability-mass function (the weak-
derivative branching step) indexes into this order, so two samplers must present
the same clocks in the same positions.

Implemented ONCE, generically, over the per-sampler accessor
[`enabling_times`](@ref): the query needs nothing a correct generalized
semi-Markov backend was not already storing. The caller turns the ages into a
selection pmf (and its θ-derivative) by rebuilding distributions at dual θ and
evaluating hazards at these ages.

Throws `ArgumentError` on a sampler whose `supports_enabled_ages` is `false`.
"""
function enabled_ages(s::SSA{K,T}, when::T) where {K,T}
    supports_enabled_ages(s) || throw(ArgumentError(
        "enabled_ages is not supported by $(typeof(s)); " *
        "supports_enabled_ages(::$(typeof(s))) is false. An exponential-only " *
        "sampler keeps no enabling-time table, so age is undefined for it."))
    ages = Tuple{K,Float64}[(key, Float64(when - te)) for (key, te) in enabling_times(s)]
    sort!(ages; by=first)
    return ages
end


"""
    retained_draw(sampler, key) -> (te=te, u=u, logu=logu)

The enabling time `te` and the SURVIVAL-SPACE uniform `u` behind clock `key`'s
current tentative firing time, satisfying the retained-draw identity

    tentative_time == te + invlogccdf(dist, logu),   u == exp(logu).

`logu` (the log-survival coordinate) is returned alongside `u` because the
identity lives in log space: a deep tail underflows `u` to `0.0` while `logu`
stays finite and inverts back to the firing time without loss. The value is a
`NamedTuple` so callers read fields by name.

Throws `ArgumentError` on a sampler whose `supports_retained_draw` is `false`.
Only `CombinedNextReaction` retains the draw randomness this reports.
"""
function retained_draw(s::SSA, key)
    throw(ArgumentError(
        "retained_draw is not supported by $(typeof(s)); " *
        "supports_retained_draw(::$(typeof(s))) is false. Only a scheduling " *
        "backend that retains draw randomness (CombinedNextReaction) can " *
        "return a retained uniform."))
end


"""
    force_fire!(sampler, key, tstar)

Fire the CHOSEN clock `key` at the CHOSEN time `tstar`, regardless of which
clock would have won the race — the branch step of the weak-derivative
estimator. The fired clock is consumed exactly as `fire!` consumes it; every
surviving clock ends up distributed by its lifetime law conditioned on survival
past `tstar`.

No random number generator is passed. When a scheduling backend must redraw a
loser whose promised time `tstar` overran, it draws from that LOSER's own keyed
stream — which is exactly what makes the forced redraw identical across a coupled
clone pair, so a coupled θ / θ+h experiment stays coupled through a forced firing.

UNBIASEDNESS PRECONDITION: the keep-if-later repair a scheduling backend applies
is proven correct when `tstar` is the current race's own decision time, which is
how the branching estimator always calls it. Choosing `tstar` by peeking at a
survivor's stored schedule biases the "keep" branch.

Throws `ArgumentError` on a sampler whose `supports_force` is `false`.
"""
function force_fire!(s::SSA{K,T}, key::K, tstar::T) where {K,T}
    throw(ArgumentError(
        "force_fire! is not supported by $(typeof(s)); " *
        "supports_force(::$(typeof(s))) is false."))
end
