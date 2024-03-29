using Random


mutable struct CommonRandomRecorder{Sampler,T,RNG}
    sampler::Sampler
    record::Dict{T,Array{RNG,1}} # Value of RNG for each clock each instance.
    sample_index::Dict{T,Int} # Current number of times each clock seen.
    miss::Dict{T,Int} # Number of misses for each clock.
end
export CommonRandomRecorder
export reset!

"""
Common random variates, also called common random numbers (CRN),
are a technique for reducing variance when repeating runs
of simulations with different parameters. The idea is to record the randomness of
one simulation and replay the same choices in subsequent runs. This particular
implementation does this by saving the state of the random number generator every
time it's used by a sampler.

The Xoshiro sampler has a relatively small state (32 bytes), which is saved
every time the sampler uses random numbers. This CRN recorder saves data in memory,
but we could save that to a memory-mapped file so that the operating system
will optimize transfer of that memory to disk.

What happens when replays of simulation runs use more draws than the first,
recorded simulation? Those simulations draw from a fresh random number
generator. This is not an exact approach.

# Example

The goal is to run the simulation with ten different parameter sets and measure
how much different parameters change the mean of some quantity determined by
the trajectories.

```julia
using Random: Xoshiro
using Fleck
example_clock = (3, 7)  # We will use clock IDs that are a tuple of 2 integers.
sampler = FirstToFire{typeof(example_clock)}()
crn_sampler = CommonRandomRecorder(sampler, typeof(example_clock), Xoshiro)
for trial_idx in 1:100
    run_simulation(model, crn_sampler)
    reset!(crn_sampler)
end
for param_idx in 1:10
    each_model = modify_model!(model, param_idx)
    run_simulation(each_model, crn_sampler)
    reset!(crn_sampler)
end
```

"""
function CommonRandomRecorder(sampler::Sampler, clock_type::Type, rng_type::Type) where {Sampler}
    storage = Dict{clock_type,Array{rng_type,1}}()
    index = Dict{clock_type,Int64}()
    miss = Dict{clock_type,Int}()
    CommonRandomRecorder{Sampler,clock_type,rng_type}(sampler, storage, index, miss)
end


"""
    reset!(recorder::CommonRandomRecorder)

The common random recorder records the state of the random number generator
for each clock, but the same clock can be enabled multiple times in one
simulation, so it records the generator state for each (clock, index of
the enabling of that clock). The `reset!` function says we are starting a new
simulation run, so all clocks haven't been seen.
"""
function reset!(recorder::CommonRandomRecorder{S,K,R}) where {S,K,R}
    reset!(recorder.sampler)
    recorder.sample_index = Dict{K,Int64}()
    recorder.miss = Dict{K,Int}()
end

"""
    misscount(recorder::CommonRandomRecorder)

The common random recorder watches a simulation and replays the states of the
random number generator on subsequent runs. This counts the number of times
during the most recent run that a clock event happened that could not be
replayed.
"""
misscount(recorder::CommonRandomRecorder) = sum(values(recorder.miss))

"""
    misses(recorder::CommonRandomRecorder)

This iterates over pairs of misses in the common random recorder during
the most recent simulation run, where the start of a simulation run
was marked by calling `reset!`.
"""
misses(recorder::CommonRandomRecorder) = pairs(recorder.miss)

# If you weren't using clock IDs and were instead using a direct method
# it would be better to instrument the next function and record every draw
# in a vector. So it makes sense to have different CRN implementations
# for different purposes.
function next(cr::CommonRandomRecorder{Sampler}, when, rng::AbstractRNG) where {Sampler}
    return next(cr.sampler, when, rng)
end


function enable!(
    cr::CommonRandomRecorder{Sampler}, clock::T, distribution::UnivariateDistribution,
    te, when, rng::AbstractRNG) where {Sampler, T}

    clock_seen_cnt = 1 + get(cr.sample_index, clock, 0)
    samples = get(cr.record, clock, nothing)
    if samples !== nothing && clock_seen_cnt <= length(samples)
        saved_rng = copy(samples[clock_seen_cnt])
        enable!(cr.sampler, clock, distribution, te, when, saved_rng)
        # Only increment the seen count if the rng was used because Next Reaction avoids RNG.
        if saved_rng != samples[clock_seen_cnt]
            cr.sample_index[clock] = clock_seen_cnt
        end  # else don't increment.
    else
        rng_save = copy(rng)
        enable!(cr.sampler, clock, distribution, te, when, rng)
        if rng != rng_save
            if samples !== nothing
                push!(samples, rng_save)
            else
                cr.record[clock] = [rng_save]
            end
            cr.sample_index[clock] = clock_seen_cnt
        end
        cr.miss[clock] = get(cr.miss, clock, 0) + 1
    end
end


function disable!(cr::CommonRandomRecorder{Sampler}, clock::T, when) where {Sampler, T}
    disable!(cr.sampler, clock, when)
end
