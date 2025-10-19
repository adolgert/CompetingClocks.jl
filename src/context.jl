export SamplerBuilder, available_samplers, add_group!, build_sampler
export SamplingContext, enable!, fire!, isenabled, freeze!

has_steploglikelihood(::Type) = false
has_steploglikelihood(::Type{<:CombinedNextReaction}) = true
has_steploglikelihood(::Type{<:DirectCall}) = true
has_steploglikelihood(::Type{<:EnabledWatcher}) = true


mutable struct SamplerBuilderGroup
    name::Symbol
    selector::Union{Function,Nothing}
    sampler_spec::Tuple{Symbol}
    instance::SSA
    SamplerBuilderGroup(name::Symbol, selector, sampler_spec) = new(name, selector, sampler_spec)
end



struct SamplerBuilder{K,T}
    clock_type::Type{K}
    time_type::Type{T}
    step_likelihood::Bool
    trajectory_likelihood::Bool
    debug::Bool
    recording::Bool
    common_random::Bool
    group::Vector{SamplerBuilderGroup}
    samplers::Dict{Tuple{Symbol,Vararg{Symbol}},Function}
    start_time::T
end

"""
    SamplerBuilder(::Type{K}, ::Type{T};
        step_likelihood=false,
        trajectory_likelihood=false,
        debug=false,
        recording=false,
        common_random=false,
        sampler_spec=:none)

A SamplerBuilder is responsible for recording a user's requirements and building
an initial sampler.

 * `K` and `T` are the clock type and time type.
 * `step_likelihood` - whether you will call `steploglikelihood` before each `fire!`
 * `trajectory_likelihood` - whether you will call `trajectoryloglikelihood`
    at the end of a simulation run.
 * `debug` - Print log messages at the debug level.
 * `recording` - Store every enable and disable for later examination.
 * `common_random` - Use common random numbers during sampling.
 * `sampler_spec` - If you want a single, particular sampler, put its Symbol name here.

# Example

```julia
builder = SamplerBuilder(Tuple,Float64)
add_group!(builder, :sparky => (x,d) -> x[1] == :recover, sampler_spec=(:nextreaction,))
add_group!(builder, :forthright=>(x,d) -> x[1] == :infect)
context = SamplingContext(builder, rng)
```
"""
function SamplerBuilder(::Type{K}, ::Type{T};
    step_likelihood=false,
    trajectory_likelihood=false,
    debug=false,
    recording=false,
    common_random=false,
    sampler_spec::Union{Symbol,Tuple{Symbol}}=(:none,),    # Ask for specific sampler.
    start_time::T=zero(T),
) where {K,T}
    group = SamplerBuilderGroup[]
    avail = make_builder_dict()
    builder = SamplerBuilder(
        K, T, step_likelihood, trajectory_likelihood, debug, recording,
        common_random, group, avail, start_time
    )
    if sampler_spec != (:none,)
        add_group!(builder, :all => (x, d) -> true; sampler_spec=sampler_spec)
    end
    return builder
end

available_samplers(builder::SamplerBuilder) = keys(builder.samplers)


"""
The `selector` defines the group of clocks that go to this sampler using
an inclusion rule, so it's a function from a clock key and distribution to a Bool.
"""
function add_group!(
    builder::SamplerBuilder,
    selector::Union{Pair,Nothing}=nothing;      # Which clocks use this sampler.
    sampler_spec::Union{Symbol,Tuple{Symbol}}=(:any,),    # Ask for specific sampler.
)
    sampler_spec = sampler_spec isa Symbol ? (sampler_spec,) : sampler_spec
    sampler_spec = sampler_spec == (:any,) ? (:firsttofire,) : sampler_spec
    if sampler_spec ∉ keys(builder.samplers)
        error("Looking for a sampler in this list: $(keys(builder.samplers))")
    end
    if length(builder.group) >= 1 && (builder.group[1].selector === nothing || selector === nothing)
        error("Need a selector on all samplers if there is more than one sampler.")
    end
    name = selector isa Pair ? selector.first : :all
    selector_func = selector isa Pair ? selector.second : selector
    push!(builder.group, SamplerBuilderGroup(name, selector_func, sampler_spec))
    return nothing
end


function make_builder_dict()
    return Dict([
        (:nextreaction,) => (K, T) -> CombinedNextReaction{K,T}(),
        (:direct,) => (K, T) -> DirectCallExplicit(K, T, KeyedRemovalPrefixSearch, BinaryTreePrefixSearch),
        (:direct, :remove, :tree) => (K, T) -> DirectCallExplicit(K, T, KeyedRemovalPrefixSearch, BinaryTreePrefixSearch),
        (:direct, :keep, :tree) => (K, T) -> DirectCallExplicit(K, T, KeyedKeepPrefixSearch, BinaryTreePrefixSearch),
        (:direct, :remove, :array) => (K, T) -> DirectCallExplicit(K, T, KeyedRemovalPrefixSearch, CumSumPrefixSearch),
        (:direct, :keep, :array) => (K, T) -> DirectCallExplicit(K, T, KeyedKeepPrefixSearch, CumSumPrefixSearch),
        (:firstreaction,) => (K, T) -> FirstReaction{K,T}(),
        (:firsttofire,) => (K, T) -> FirstToFire{K,T}(),
        (:petri,) => (K, T) -> Petri{K,T}(),
    ])
end

"""
See sampler.jl for the MultiSampler to understand how we're making a chooser.
We would like to implement this:
"Compiling Pattern Matching to Good Decision Trees" - Luc Maranget (2008).
It looks like MLStyle.jl has a good example of how to do it.
"""
struct FromInclusion{K} <: SamplerChoice{Symbol,K}
    matcher::Function
end

function CompetingClocks.choose_sampler(
    chooser::FromInclusion, clock, distribution::UnivariateDistribution)
    return chooser.matcher(clock, distribution)
end

function build_sampler(builder::SamplerBuilder)
    isempty(builder.group) && error("Need to add_group! on the builder.")
    K = builder.clock_type
    T = builder.time_type
    if length(builder.group) == 1
        sampler = builder.samplers[builder.group[1].sampler_spec](K, T)
        matcher = nothing
    else
        # All Direct methods get combined into one.
        directs = filter(x -> x.sampler_spec[1] == :direct, builder.group)
        if length(directs) > 1
            error("implement soon for multiple directs")
        elseif length(directs) == 1
            direct = directs[1]
            direct.instance = builder.samplers[direct.sampler_spec](K, T)
        else
            direct = nothing
        end
        competes = filter(x -> x.sampler_spec[1] != :direct, builder.group)
        for compete in competes
            compete.instance = builder.samplers[compete.sampler_spec](K, T)
        end
        # Any direct method gets added to the others for combination.
        if !isnothing(direct)
            push!(competes, direct)
        end
        inclusion = Dict(samp.name => samp.selector for samp in competes)
        matcher = FromInclusion{K}(make_key_classifier(inclusion))
        sampler = MultiSampler{Symbol,K,T}(matcher)
        for samp in competes
            sampler[samp.name] = samp.instance
        end
    end
    return sampler
end


function make_key_classifier(inclusion_criteria::Dict{Symbol,Function})
    # Convert to tuple for better type stability
    criteria_tuple = Tuple(inclusion_criteria)

    return function (key, dist)
        for (symbol, criterion_func) in criteria_tuple
            if criterion_func(key, dist)
                return symbol
            end
        end
        error("No matching set found for key: $key")
    end
end


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
mutable struct SamplingContext{K,T,Sampler<:SSA{K,T},RNG,Like,CRN,Dbg}
    sampler::Sampler # The actual sampler
    rng::RNG
    likelihood::Like
    crn::CRN
    debug::Dbg      # Union{Nothing, TrackingState{K,T}}
    split_weight::Float64
    time::T
    fixed_start::T
end


function SamplingContext(builder::SamplerBuilder, rng::R) where {R<:AbstractRNG}
    K = builder.clock_type
    T = builder.time_type
    sampler = build_sampler(builder)
    if builder.trajectory_likelihood
        likelihood = TrajectoryWatcher{K,T}()
    else
        likelihood = nothing
    end
    if builder.step_likelihood && !has_steploglikelihood(typeof(sampler)) && isnothing(likelihood)
        likelihood = TrackWatcher{K,T}()
    end
    if builder.common_random
        crn = CommonRandom{K,R}()
    else
        crn = nothing
    end
    if builder.debug || builder.recording
        debug = DebugWatcher{K,T}(log=builder.debug)
    else
        debug = nothing
    end
    SamplingContext{K,T,typeof(sampler),R,typeof(likelihood),typeof(crn),typeof(debug)}(
        sampler, rng, likelihood, crn, debug, 1.0, builder.start_time, builder.start_time
    )
end


Base.time(ctx::SamplingContext) = ctx.time

function freeze!(ctx::SamplingContext)
    if ctx.crn !== nothing
        freeze!(ctx.crn)
        ctx.time = ctx.fixed_start
    else
        error("Ask for common random numbers when creating the builder.")
    end
end


function enable!(ctx::SamplingContext{K,T}, clock::K, dist, relative_te::T) where {K,T}
    when = ctx.time
    te = when + relative_te
    ctx.likelihood !== nothing && enable!(ctx.likelihood, clock, dist, te, when, ctx.rng)
    if ctx.crn !== nothing
        with_common_rng(ctx.crn, clock, ctx.rng) do wrapped_rng
            enable!(ctx.sampler, clock, dist, te, when, wrapped_rng)
        end
    else
        enable!(ctx.sampler, clock, dist, te, when, ctx.rng)
    end
    ctx.debug !== nothing && enable!(ctx.debug, clock, dist, te, when, ctx.rng)
end

function enable!(ctx::SamplingContext{K,T}, clock::K, dist) where {K,T}
    enable!(ctx, clock, dist, zero(T))
end

function disable!(ctx::SamplingContext{K,T}, clock::K) where {K,T}
    when = ctx.time
    ctx.likelihood !== nothing && disable!(ctx.likelihood, clock, when)
    disable!(ctx.sampler, clock, when)
    ctx.debug !== nothing && disable!(ctx.debug, clock, when)
end


function next(ctx::SamplingContext{K,T}) where {K,T}
    next(ctx.sampler, ctx.time, ctx.rng)
end


function fire!(ctx::SamplingContext{K,T}, clock::K, when::T) where {K,T}
    ctx.likelihood !== nothing && fire!(ctx.likelihood, clock, when)
    fire!(ctx.sampler, clock, when)
    ctx.debug !== nothing && fire!(ctx.debug, clock, when)
    ctx.time = when
end


function reset!(ctx::SamplingContext{K,T}) where {K,T}
    ctx.likelihood !== nothing && reset!(ctx.likelihood)
    reset!(ctx.sampler)
    ctx.debug !== nothing && reset!(ctx.debug)
    ctx.time = ctx.fixed_start
end


function Base.copy!(dst::SamplingContext{K,T}, src::SamplingContext{K,T}) where {K,T}
    copy!(dst.sampler, src.sampler)
    copy!(dst.rng, src.rng) # Will need to initialize this separately.
    src.likelihood !== nothing && copy!(dst.likelihood, src.likelihood)
    src.debug !== nothing && copy!(dst.debug, src.debug)
    dst.split_weight = src.split_weight
    return dst
end

function enabled(ctx::SamplingContext)
    ctx.likelihood !== nothing && return enabled(ctx.likelihood)
    return enabled(ctx.sampler)
end


function Base.length(ctx::SamplingContext)
    ctx.likelihood !== nothing && return length(ctx.likelihood)
    return length(ctx.sampler)
end


function isenabled(ctx::SamplingContext{K}, clock::K) where {K}
    ctx.likelihood !== nothing && return isenabled(ctx.likelihood, clock)
    return isenabled(ctx.sampler, clock)
end

Base.keytype(ctx::SamplingContext{K}) where {K} = K
timetype(ctx::SamplingContext{K,T}) where {K,T} = T


function steploglikelihood(dc::SamplingContext, when, which)
    return steploglikelihood(ctx.sampler, ctx.when, when, which)
end


function trajectoryloglikelihood(dc::SamplingContext)
    isnothing(ctx.likelihood) && error(
        "Must enabled trajectory likelihood when creating sampler"
    )
    return trajectoryloglikelihood(ctx.likelihood)
end
