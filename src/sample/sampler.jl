
export MultiSampler


abstract type SamplerChoice{SamplerKey,Key} end

function choose_sampler(
    chooser::SamplerChoice{SamplerKey,Key}, clock::Key, distribution::UnivariateDistribution
)::SamplerKey where {SamplerKey,Key}
    throw(MissingException("No sampler choice given to the MultiSampler"))
end

export SamplerChoice
export choose_sampler

"""
    MultiSampler{SamplerKey,Key,Time}(which_sampler::Chooser) <: SSA{Key,Time}

A sampler returns the soonest event, so we can make a hierarchical sampler
that returns the soonest event of the samplers it contains. This is useful because
the performance of a sampler depends on the type of the event. For instance,
some simulations have a few fast events and a lot of slow ones, so it helps
to split them into separate data structures.

The `SamplerKey` is the type of an identifier for the samplers that this
`MultiSampler` contains. The `which_sampler` argument is a strategy object
that decides which event is sampled by which contained sampler. There is
an example of this below.

Once a clock is first enabled, it will always go to the same sampler.
This sampler remembers the associations, which could increase memory for
simulations with semi-infinite clocks.

# Examples

Let's make one sampler for exponential distributions, one for a few
clocks we know will be fast and one for slower clocks.
We can name them with symbols. The trick is that we need to direct each kind
of distribution to the correct sampler. Use a Float64 for time and each clock
can be identified with an Int64.
```
using CompetingClocks
using Distributions: Exponential, UnivariateDistribution

struct ByDistribution <: SamplerChoice{Int64,Symbol} end

function CompetingClocks.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::Exponential
    )::Symbol
    return :direct
end
function CompetingClocks.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::UnivariateDistribution
    )::Symbol
    if clock < 100
        return :fast
    else
        return :slow
    end
end
sampler = MultiSampler{Symbol,Int64,Float64}(ByDistribution())
sampler[:direct] = OptimizedDirect{Int64,Float64}()
sampler[:fast] = FirstToFire{Int64,Float64}()
sampler[:slow] = FirstToFire{Int64,Float64}()
```

"""
mutable struct MultiSampler{SamplerKey,Key,Time,Chooser} <: SSA{Key,Time}
    propagator::Dict{SamplerKey,SSA{Key,Time}}
    chooser::Chooser
    chosen::Dict{Key,SamplerKey}
end


function MultiSampler{SamplerKey,Key,Time}(
    which_sampler::Chooser
) where {SamplerKey,Key,Time,Chooser<:SamplerChoice{SamplerKey,Key}}

    MultiSampler{SamplerKey,Key,Time,Chooser}(
        Dict{SamplerKey,SSA{Key,Time}}(),
        which_sampler,
        Dict{Key,SamplerKey}()
    )
end


function reset!(sampler::MultiSampler)
    for clear_sampler in values(sampler.propagator)
        reset!(clear_sampler)
    end
    empty!(sampler.chosen)
end


function Base.copy!(
    dst::MultiSampler{SamplerKey,Key,Time,Chooser},
    src::MultiSampler{SamplerKey,Key,Time,Chooser}
) where {SamplerKey,Key,Time,Chooser}

    copy!(dst.propagator, src.propagator)
    dst.chooser = src.chooser
    copy!(dst.chosen, src.chosen)
    dst
end


function Base.setindex!(
    sampler::MultiSampler{SamplerKey,Key,Time}, algorithm::SSA{Key,Time}, sampler_key::SamplerKey
) where {SamplerKey,Key,Time}
    sampler.propagator[sampler_key] = algorithm
end


function next(
    sampler::MultiSampler{SamplerKey,Key,Time},
    when::Time,
    rng::AbstractRNG
) where {SamplerKey,Key,Time}

    least_when::Time = typemax(Time)
    least_transition::Union{Nothing,Key} = nothing
    for propagator in values(sampler.propagator)
        when, transition = next(propagator, when, rng)
        if when < least_when
            least_when = when
            least_transition = transition
        end
    end
    return (least_when, least_transition)
end


function enable!(
    sampler::MultiSampler{SamplerKey,Key,Time},
    clock::Key,
    distribution::UnivariateDistribution,
    te::Time,
    when::Time,
    rng::AbstractRNG
) where {SamplerKey,Key,Time}
    @debug "Enabling the MultiSampler"
    this_clock_sampler = choose_sampler(sampler.chooser, clock, distribution)
    sampler.chosen[clock] = this_clock_sampler
    propagator = sampler.propagator[this_clock_sampler]
    enable!(propagator, clock, distribution, te, when, rng)
end


function disable!(
    sampler::MultiSampler{SamplerKey,Key,Time}, clock::Key, when::Time
) where {SamplerKey,Key,Time}
    disable!(sampler.propagator[sampler.chosen[clock]], clock, when)
end


function Base.getindex(sampler::MultiSampler, clock)
    return getindex(sampler.chosen[clock], clock)
end


function Base.keys(sampler::MultiSampler)
    return union([keys(propagator) for propagator in values(sampler.propagator)]...)
end


function Base.length(sampler::MultiSampler)
    return sum([length(propagator) for propagator in values(sampler.propagator)])
end

function Base.haskey(sampler::MultiSampler{SamplerKey,Key,Time,Chooser}, clock) where {SamplerKey,Key,Time,Chooser}
    for propagator in values(sampler.propagator)
        if haskey(propagator, clock)
            return true
        else
            continue
        end
    end
    return false
end
