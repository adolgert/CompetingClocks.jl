using Random: AbstractRNG
using Distributions: UnivariateDistribution

export enable!, disable!, next


function enable!(
    sampler,
    clock,
    distribution::UnivariateDistribution,
    te::Float64, # enabling time
    when::Float64, # current simulation time
    rng::AbstractRNG
    )
end

function disable!(sampler, clock, when::Float64) end

function next(sampler, when::Float64, rng::AbstractRNG) end
