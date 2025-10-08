"""
The SamplingContext is responsible for doing a perfect forwarding to multiple
components that implement the desired features in a sampler.

It uses a Context pattern where each type parameter is either present or
Nothing, and the compiler knows to optimize-away the Nothing.
We keep the internal logic simple. If something is present, call it.

 - memory
 - delayed transitions
 - hierarchical samplers
"""
struct SamplingContext{K,T,Sampler<:SSA{K,T},RNG,Rec,Dbg}
    sampler::Sampler # The actual sampler
    rng::RNG
    recording::Rec  # Union{Nothing, RecordingState{K,RNG}}
    debug::Dbg      # Union{Nothing, TrackingState{K,T}}
    split_weight::Float64
end


step_likelihood_api(::Type) = false
track_likelihood_api(::Type) = false
enabled_api(::Type) = false

#=
M = maybe soon. X = yes. ' ' = no.
FTF-S is a first-to-fire version that stores more information.

        STEP TRACK ENABLED
NR      M    M     X
Direct  X    M     X
FR      X          X
FR-T    X    X     X
FTF
FTF-S   X          X
FTF-ST  X    X     X
Petri   X          X
Petri-T X    X     X
=#

function make_sampler(
    ;
)
end


function enable!(ctx::SamplingContext{K,T}, clock::K, dist, te::T, when::T) where {K,T}
    enable!(ctx.sampler, clock, dist, te, when, ctx.rng)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


function disable!(ctx::SamplingContext{K,T}, clock::K, when::T) where {K,T}
    disable!(ctx.sampler, clock, when)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


function next(ctx::SamplingContext{K,T}, when::T) where {K,T}
    next(ctx.sampler, when, ctx.rng)
end


function fire!(ctx::SamplingContext{K,T}, clock::K, dist, te::T, when::T) where {K,T}
    fire!(ctx.sampler, clock, dist, te, when, ctx.rng)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_fire!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_fire!(ctx.debug, clock, dist, te, when)
end

function reset!(ctx::SamplingContext{K,T}) where {K,T}
    reset!(ctx.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


function Base.copy!(src::SamplingContext{K,T}, dst::SamplingContext{K,T}) where {K,T}
    copy!(src.sampler, dst.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    src.recording !== nothing && record_copy!(src.recording, clock, src.rng)
    src.tracking !== nothing && track_copy!(src.tracking, clock, dist, te, when)
end


function Base.getindex(ctx::SamplingContext{K}, clock::K) where {K}
    # Core sampler logic first
    getindex(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


function Base.keys(ctx::SamplingContext)
    # Core sampler logic first
    keys(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_keys(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_keys(ctx.debug, clock, dist, te, when)
end


function Base.length(ctx::SamplingContext)
    # Core sampler logic first
    length(ctx.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


function Base.haskey(ctx::SamplingContext{K}, clock::K) where {K}
    # Core sampler logic first
    haskey(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.debug !== nothing && track_enable!(ctx.debug, clock, dist, te, when)
end


# Get key type (already defined, not overridden)
Base.keytype(ctx::SamplingContext{K}) where {K} = K

# Get time type
timetype(ctx::SamplingContext{K,T}) where {K,T} = T

function steploglikelihood(dc::SamplingContext, now, when, which)
    return steploglikelihood(ctx.sampler, now, when, which)
end


function trajectoryloglikelihood(dc::SamplingContext)
    return trajectoryloglikelihood(ctx.sampler)
end
