using Random


"""
Common random variates are a technique for reducing variance when repeating runs
of simulations with different parameters. The idea is to record the randomness of
one simulation and replay the same choices in subsequent runs. This particular
implementation does this by saving the state of the random number generator every
time it's used by a sampler.

```julia
using Random: Xoshiro
example_clock = (3, 7)
sampler = FirstToFire{typeof(example_clock)}()
crn_sampler = CommonRandomRecorder(sampler, typeof(example_clock), Xoshiro)
run_simulation(model, crn_sampler)
for param_idx in 1:10
    each_model = modify_model(model, param_idx)
    replay_sampler = replay(crn_sampler)
    run_simulation(each_model, replay_sampler)
end
```

The Xoshiro sampler has a relatively small state (32 bytes), which is saved
every time the sampler uses random numbers. This version saves data in memory.
"""
struct CommonRandomRecorder{Sampler,T,RNG}
    sampler::Sampler
    record::Dict{T,Array{RNG,1}}
end
export CommonRandomRecorder

function CommonRandomRecorder(sampler::Sampler, clock_type::Type, rng_type::Type) where {Sampler}
    storage = Dict{clock_type,Array{rng_type,1}}()
    CommonRandomRecorder{Sampler,clock_type,rng_type}(sampler, storage)
end


# If you weren't using clock IDs and were instead using a direct method
# it would be better to instrument the next function and record every draw
# in a vector. So it makes sense to have different CRN implementations
# for different purposes.
function next(cr::CommonRandomRecorder{Sampler}, when::Float64, rng::AbstractRNG) where {Sampler}
    return next(cr.sampler, when, rng)
end


function enable!(
    cr::CommonRandomRecorder{Sampler}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {Sampler, T}

    rng_save = copy(rng)
    enable!(cr.sampler, clock, distribution, te, when, rng)
    if rng != rng_save
        if clock ∈ keys(cr.record)
            push!(cr.record[clock], rng_save)
        else
            cr.record[clock] = [rng_save]
        end
    end
end


function disable!(cr::CommonRandomRecorder{Sampler}, clock::T, when::Float64) where {Sampler, T}
    disable!(cr.sampler, clock, when)
end


struct CommonRandomReplay{CommonRandomRecorder,T}
    recorder::CommonRandomRecorder
    sample_index::Dict{T,Int}
    miss::Dict{T,Int}
end


function replay(recorder::CommonRandomRecorder{Sampler,T,RNG}) where {Sampler,T,RNG}
    index = Dict(clock => 1 for clock in keys(recorder.record))
    CommonRandomReplay{CommonRandomRecorder{Sampler,T,RNG},T}(recorder, index, Dict{T,Int}())
end


function with_generator(fn::Function, crr::CommonRandomReplay, clock, rng)
    clock_index = get(crr.sample_index, clock, 0)
    if clock_index > 0
        samples = crr.recorder.record[clock]
        if clock_index <= length(samples)
            use_rng = copy(samples[clock_index])
            fn(use_rng)
            if use_rng != samples[clock_index]
                crr.sample_index[clock] += 1
            end
            return
        end
    end
    fn(rng)
    crr.miss[clock] = get(crr.miss, clock, 0) + 1
end


function next(crr::CommonRandomReplay{Sampler}, when::Float64, rng::AbstractRNG) where {Sampler}
    return next(crr.sampler, when, rng)
end


function enable!(
    crr::CommonRandomReplay{Sampler}, clock::T, distribution::UnivariateDistribution,
    te::Float64, when::Float64, rng::AbstractRNG) where {Sampler, T}

    with_generator(crr, clock, rng) do use_rng
        enable!(crr.sampler, clock, distribution, te, when, use_rng)
    end
end


function disable!(cr::CommonRandomReplay{Sampler}, clock::T, when::Float64) where {Sampler, T}
    disable!(crr.sampler, clock, when)
end
