
struct SamplingContext{K,T,Sampler<:SSA{K,T},RNG,Rec,Trk,Lk}
    sampler::Sampler # The actual sampler
    rng::RNG
    recording::Rec  # Union{Nothing, RecordingState{K,RNG}}
    tracking::Trk   # Union{Nothing, TrackingState{K,T}}
    likelihood::Lk  # Union{Nothing, LikelihoodState{K,T}}
    split_weight::Float64
end


function enable!(ctx::SamplingContext{K,T}, clock::K, dist, te::T, when::T) where {K,T}
    # Core sampler logic first
    enable!(ctx.sampler, clock, dist, te, when, ctx.rng)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function disable!(ctx::SamplingContext{K,T}, clock::K, when::T) where {K,T}
    # Core sampler logic first
    disable!(ctx.sampler, clock, when)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function next(ctx::SamplingContext{K,T}, when::T) where {K,T}
    # Core sampler logic first
    next(ctx.sampler, when, ctx.rng)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function reset!(ctx::SamplingContext{K,T}) where {K,T}
    # Core sampler logic first
    reset!(ctx.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function Base.copy!(src::SamplingContext{K,T}, dst::SamplingContext{K,T}) where {K,T}
    # Core sampler logic first
    copy!(src.sampler, dst.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    src.recording !== nothing && record_enable!(src.recording, clock, src.rng)
    src.tracking !== nothing && track_enable!(src.tracking, clock, dist, te, when)
    src.likelihood !== nothing && likelihood_enable!(src.likelihood, clock, dist, te, when)
end


function Base.getindex(ctx::SamplingContext{K}, clock::K) where {K}
    # Core sampler logic first
    getindex(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function Base.keys(ctx::SamplingContext)
    # Core sampler logic first
    keys(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function Base.length(ctx::SamplingContext)
    # Core sampler logic first
    length(ctx.sampler)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


function Base.haskey(ctx::SamplingContext{K}, clock::K) where {K}
    # Core sampler logic first
    haskey(ctx.sampler, clock)

    # Feature hooks (compiler eliminates branches when types are Nothing)
    ctx.recording !== nothing && record_enable!(ctx.recording, clock, ctx.rng)
    ctx.tracking !== nothing && track_enable!(ctx.tracking, clock, dist, te, when)
    ctx.likelihood !== nothing && likelihood_enable!(ctx.likelihood, clock, dist, te, when)
end


# Get key type (already defined, not overridden)
Base.keytype(ctx::SamplingContext{K}) where {K} = K

# Get time type
timetype(ctx::SamplingContext{K,T}) where {K,T} = T

function steploglikelihood(dc::SamplingContext, now, when, which)
end
