using Random

abstract type CommonRandom{S,K,R} end

struct CommonRandomRecorder{Sampler,K,RNG} <: CommonRandom{Sampler,K,RNG}
    sampler::Sampler
    record::Dict{K,Array{RNG,1}} # Value of RNG for each clock each instance.
    sample_index::Dict{K,Int} # Current number of times each clock seen.
    miss::Dict{K,Int} # Number of misses for each clock.
end
export CommonRandomRecorder, FrozenCommonRandomRecorder
export reset!, misscount, misses, freeze

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
using CompetingClocks
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
    empty!(recorder.sample_index)
    empty!(recorder.miss)
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
function next(cr::CommonRandomRecorder, when, rng::AbstractRNG)
    return next(cr.sampler, when, rng)
end


function enable!(
    cr::CommonRandomRecorder, clock, distribution::UnivariateDistribution,
    te, when, rng::AbstractRNG)

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


function disable!(cr::CommonRandomRecorder, clock, when)
    disable!(cr.sampler, clock, when)
end


struct FrozenCommonRandomRecorder{S,K,R} <: CommonRandom{S,K,R}
    cr::CommonRandomRecorder{S,K,R}
    miss::Dict{K,Int}
end


"""
    freeze(cr::CommonRandomRecorder)::FrozenCommonRandomRecorder

The [CommonRandomRecorder](@ref) records every time it sees a clock request
random number generation. It continues to do that every time it runs, which
is a problem if you run simulations for comparison on multiple threads.
If you want to use CRN and to use multiple threads for subsequent simulation
runs, then first run the simulation a bunch of times on one thread. Then
freeze the simulation, and then the frozen version will stop remembering
new threads.

There is one part of the frozen recorder that will be mutable because it's
useful for debugging, the record of missed clocks. Freeze a recorder for
each thread, and each thread will track its own misses. They will all work
from the same copy of the recorded random number generator states.
"""
freeze(cr::CommonRandomRecorder{S,K,R}) where {S,K,R} = FrozenCommonRandomRecorder(cr, Dict{K,Int}())
reset!(fcr::FrozenCommonRandomRecorder{S,K,R}) where {S,K,R} = (reset!(fcr.cr); empty!(fcr.miss); nothing)
misscount(fcr::FrozenCommonRandomRecorder) = misscount(cr.miss)
misses(fcr::FrozenCommonRandomRecorder) = misses(fcr.miss)
next(fcr::FrozenCommonRandomRecorder, when, rng::AbstractRNG) = next(fcr.cr, when, rng)
disable!(fcr::FrozenCommonRandomRecorder, clock, when) = disable!(fcr.cr, clock, when)

function enable!(
    fcr::FrozenCommonRandomRecorder, clock, distribution::UnivariateDistribution,
    te, when, rng::AbstractRNG)

    cr = fcr.cr
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
        enable!(cr.sampler, clock, distribution, te, when, rng)
        fcr.miss[clock] = get(fcr.miss, clock, 0) + 1
    end
end
