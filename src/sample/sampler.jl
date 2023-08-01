
export SingleSampler, MultiSampler
export enable!, disable!, sample!

"""
    SSA{Key,Time}

This abstract type represents a stochastic simulation algorithm (SSA). It is
parametrized by the clock ID, or key, and the type used for the time, which
is typically a Float64.
"""
abstract type SSA{Key,Time} end


"""
    SingleSampler{SSA,Time}(propagator::SSA)

This makes a sampler from a single stochastic simulation algorithm. It combines
the core algorithm with the rest of the state of the system, which is just
the time.
"""
mutable struct SingleSampler{Algorithm,Time}
    propagator::Algorithm
    when::Time
end


function SingleSampler(propagator::SSA{Key,Time}) where {Key,Time}
    SingleSampler{SSA{Key,Time},Time}(propagator, zero(Time))
end


function sample!(sampler::SingleSampler, rng::AbstractRNG)
    when, transition = next(sampler.propagator, sampler.when, rng)
    if transition !== nothing
        sampler.when = when
        disable!(sampler.propagator, transition, sampler.when)
    end
    return (when, transition)
end


function enable!(
    sampler::SingleSampler, clock, distribution::UnivariateDistribution, te, rng::AbstractRNG)
    enable!(sampler.propagator, clock, distribution, te, sampler.when, rng)
end


function disable!(sampler::SingleSampler, clock)
    disable!(sampler.propagator, clock, sampler.when)
end


abstract type SamplerChoice{Key,SamplerKey} end

function choose_sampler(
    chooser::SamplerChoice{Key,SamplerKey}, clock::Key, distribution::UnivariateDistribution
    )::SamplerKey where {Key,SamplerKey}
    throw(MissingException("No sampler choice given to the MultiSampler"))
end

export SamplerChoice
export choose_sampler

"""
    MultiSampler{SamplerKey,Key,Time}(which_sampler::Function)

This makes a sampler that uses multiple stochastic sampling algorithms (SSA) to
determine the next transition to fire. It returns the soonest transition of all
of the algorithms. The `which_sampler` function looks at the clock ID, or key,
and chooses which sampler should sample this clock. Add algorithms to this
sampler like you would add them to a dictionary.

# Examples

Let's make one sampler for exponential distributions and one for the rest.
We can name them with symbols. The trick is that we need to direct each kind
of distribution to the correct sampler. Use a Float64 for time and each clock
can be identified with an Int64.
```
using Fleck
using Distributions: Exponential, UnivariateDistribution

struct ByDistribution <: SamplerChoice{Int64,Symbol} end

function Fleck.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::Exponential
    )::Symbol
    return :direct
end
function Fleck.choose_sampler(
    chooser::ByDistribution, clock::Int64, distribution::UnivariateDistribution
    )::Symbol
    return :otherone
end
sampler = MultiSampler{Symbol,Int64,Float64}(ByDistribution())
sampler[:direct] = OptimizedDirect{Int64,Float64}()
sampler[:otherone] = FirstToFire{Int64,Float64}()
```
Why don't we choose samplers by passing a function into the `MultiSampler`?
This is an effort to ensure type safety, because if you're going to the trouble
of using a hierarchical sampler, you clearly care about speed.
"""
mutable struct MultiSampler{SamplerKey,Key,Time,Chooser}
    propagator::Dict{SamplerKey,SSA{Key,Time}}
    when::Time
    chooser::Chooser
    chosen::Dict{Key,SamplerKey}
end


function MultiSampler{SamplerKey,Key,Time}(
    which_sampler::Chooser
    ) where {SamplerKey,Key,Time,Chooser <: SamplerChoice{Key,SamplerKey}}

    MultiSampler{SamplerKey,Key,Time,Chooser}(
        Dict{SamplerKey,SSA{Key,Time}}(),
        zero(Time),
        which_sampler,
        Dict{Key,SamplerKey}()
        )
end


function Base.setindex!(
    sampler::MultiSampler{SamplerKey,Key,Time}, algorithm::SSA{Key,Time}, sampler_key::SamplerKey
    ) where {SamplerKey,Key,Time}
    sampler.propagator[sampler_key] = algorithm
end


function sample!(
    sampler::MultiSampler{SamplerKey,Key,Time},
    rng::AbstractRNG
    ) where {SamplerKey,Key,Time}

    least_when::Time = typemax(Time)
    least_transition::Union{Nothing,Key} = nothing
    least_source::Union{Nothing,SamplerKey} = nothing
    for (sample_key, propagator) in sampler.propagator
        when, transition = next(propagator, sampler.when, rng)
        if when < least_when
            least_when = when
            least_transition = transition
            least_source = sample_key
        end
    end
    if least_transition !== nothing
        sampler.when = least_when
        disable!(sampler.propagator[least_source], least_transition, least_when)
    end
    return (least_when, least_transition)
end


function enable!(
    sampler::MultiSampler, clock, distribution::UnivariateDistribution, te, rng::AbstractRNG
    )
    this_clock_sampler = choose_sampler(sampler.chooser, clock, distribution)
    sampler.chosen[clock] = this_clock_sampler
    propagator = sampler.propagator[this_clock_sampler]
    enable!(propagator, clock, distribution, te, sampler.when, rng)
end


function disable!(sampler::MultiSampler, clock)
    disable!(sampler.propagator[sampler.chosen[clock]], clock, sampler.when)
end