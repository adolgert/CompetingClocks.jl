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


struct SamplerBuilderGroup
    name::Symbol
    selector::Function
    frequency_tier::Int64
    constant::Bool
    sampler::Symbol
    space::Symbol
    order::Int64
end

struct SamplerBuilder{K,T}
    clock_type::K
    time_type::T
    step_likelihood::Bool
    trajectory_likelihood::Bool
    debug::Bool
    recording::Bool
    common_random::Bool
    group::Vector{SamplerBuilderGroup}
end

function SamplerBuilder(::Type{K}, ::Type{T};
    step_likelihood=false,
    trajectory_likelihood=false,
    debug=false,
    recording=false,
    common_random=false,
) where {K,T}
    SamplerBuilder{K,T}(K, T, step_likelihood, trajectory_likelihood, debug, recording, common_random)
end

function add_group!(
    builder::SamplerBuilder,
    name::Symbol;            # User-given name for this sampler.
    selector::Function,      # Which clocks use this sampler.
    frequency_tier::Int64=1, # Higher number = more churn.
    constant::Bool=false,    # Constant hazards, ie Exponential.
    sampler::Symbol=:any,    # Ask for specific sampler.
    space::Symbol=:countable,# :countable, :finite
    order::Int64,            # Approximate number of clocks.
)
    push!(builder.group,
        SamplerBuilderGroup(
            name, selector, frequency_tier, constant, sampler
        )
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
