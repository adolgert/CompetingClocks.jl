
export SingleSampler, MultiSampler
export enable!, disable!, sample!

"""
    SSA{Key,Time}

This abstract type represents a stochastic simulation algorithm. It is
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


"""
    MultiSampler{SamplerKey,Key,Time}(which_sampler::Function)

This makes a sampler that uses multiple stochastic sampling algorithms (SSA) to
determine the next transition to fire. It returns the soonest transition of all
of the algorithms. The `which_sampler` function looks at the clock ID, or key,
and chooses which sampler should sample this clock. Add algorithms to this
sampler like you would add them to a dictionary.
"""
mutable struct MultiSampler{SamplerKey,Key,Time}
    propagator::Dict{SamplerKey,SSA{Key,Time}}
    when::Time
    which_sampler::Function
end


function MultiSampler{SamplerKey,Key,Time}(which_sampler::Function) where {SamplerKey,Key,Time}
    MultiSampler{SamplerKey,Key,Time}(
        Dict{SamplerKey,SSA{Key,Time}}(),
        zero(Time),
        which_sampler
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
    for i in eachindex(sampler.propagator)
        when, transition = next(sampler.propagator[i], sampler.when, rng)
        if when < least_when
            least_when = when
            least_transition = transition
            least_source = i
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
    propagator = sampler.propagator[sampler.which_sampler(clock)]
    enable!(propagator, clock, distribution, te, sampler.when, rng)
end


function disable!(sampler::MultiSampler, clock)
    disable!(sampler.propagator[sampler.which_sampler(clock)], clock, sampler.when)
end
