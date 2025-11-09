export SamplerBuilder, add_group!, build_sampler

has_steploglikelihood(::Type) = false
has_steploglikelihood(::Type{<:CombinedNextReaction}) = true
has_steploglikelihood(::Type{<:DirectCall}) = true
has_steploglikelihood(::Type{<:MultipleDirect}) = true
has_steploglikelihood(::Type{<:EnabledWatcher}) = true
has_pathloglikelihood(::Type) = false
has_pathloglikelihood(::Type{DirectCall}) = true
has_pathloglikelihood(::Type{MultipleDirect}) = true


mutable struct SamplerBuilderGroup
    name::Symbol
    selector::Union{Function,Nothing}
    method::Union{SamplerSpec,Nothing}
    instance::SSA
    # Constructor sets `instance` member to undefined.
    SamplerBuilderGroup(name::Symbol, selector, method) = new(name, selector, method)
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
    start_time::T
    likelihood_cnt::Int
end

"""
    SamplerBuilder(::Type{K}, ::Type{T};
        step_likelihood=false,
        trajectory_likelihood=false,
        debug=false,
        recording=false,
        common_random=false,
        method=nothing,
        start_time::T,
        likelihood_cnt::Int
        )

A SamplerBuilder is responsible for recording a user's requirements and building
an initial sampler.

 * `K` and `T` are the clock type and time type.
 * `step_likelihood` - whether you will call `steploglikelihood` before each `fire!`
 * `trajectory_likelihood` - whether you will call `pathloglikelihood`
    at the end of a simulation run.
 * `debug` - Print log messages at the debug level.
 * `recording` - Store every enable and disable for later examination.
 * `common_random` - Use common random numbers during sampling.
 * `method` - If you want a single, particular sampler, put its `SamplerSpec` here.
   It will create a group called `:all` that has this sampling method.
 * `start_time` - Sometimes a simulation shouldn't start at zero.
 * `likelihood_cnt` - The number of likelihoods to compute, corresponds to number of
   distributions in `enable!` calls. This turns on `trajectory_likelihood`.

# Example

```julia
builder = SamplerBuilder(Tuple,Float64)
add_group!(builder, :sparky => (x,d) -> x[1] == :recover, method=NextReaction())
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
    method::Union{SamplerSpec,Nothing}=nothing,    # Ask for specific sampler.
    start_time::T=zero(T),
    likelihood_cnt=1,
) where {K,T}
    group = SamplerBuilderGroup[]
    trajectory_likelihood = trajectory_likelihood || likelihood_cnt > 1
    builder = SamplerBuilder(
        K, T, step_likelihood, trajectory_likelihood, debug, recording,
        common_random, group, start_time, likelihood_cnt
    )
    if !isnothing(method)
        add_group!(builder, :all => (x, d) -> true; method=method)
    end
    return builder
end


"""
The `selector` defines the group of clocks that go to this sampler using
an inclusion rule, so it's a function from a clock key and distribution to a Bool.
"""
function add_group!(
    builder::SamplerBuilder,
    selector::Union{Pair,Nothing}=nothing;      # Which clocks use this sampler.
    method::Union{SamplerSpec,Nothing}=nothing, # Ask for specific sampler.
)
    if length(builder.group) >= 1 && (builder.group[1].selector === nothing || selector === nothing)
        error("Need a selector on all samplers if there is more than one sampler.")
    end
    name = selector isa Pair ? selector.first : :all
    selector_func = selector isa Pair ? selector.second : selector
    push!(builder.group, SamplerBuilderGroup(name, selector_func, method))
    return nothing
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

function auto_select_method(builder::SamplerBuilder)
    # Auto-select a sampler method based on builder requirements
    if builder.trajectory_likelihood
        return DirectMethod()
    elseif builder.step_likelihood
        return NextReactionMethod()
    else
        return FirstToFireMethod()
    end
end

function build_sampler(builder::SamplerBuilder)
    K = builder.clock_type
    T = builder.time_type
    if length(builder.group) == 0
        sampler = FirstToFireMethod()(K, T)
        matcher = nothing
    elseif length(builder.group) == 1
        method = isnothing(builder.group[1].method) ? auto_select_method(builder) : builder.group[1].method
        sampler = method(K, T)
        matcher = nothing
    else
        competes = builder.group
        for compete in competes
            method = isnothing(compete.method) ? auto_select_method(builder) : compete.method
            compete.instance = method(K, T)
        end
        # Any direct method gets added to the others for combination.
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
