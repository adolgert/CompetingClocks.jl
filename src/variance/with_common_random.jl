export CommonRandom, misses, misscount, reset_crn!

"""
Continuation/Functional Pattern

  Key Insight: CommonRandom provides RNG transformation as a higher-order function.

  struct CommonRandom{K,R}
      recorder::Dict{K,Vector{R}}
      mode::Symbol
  end

  function with_common_rng(
      crc::CommonRandom{K,R}, clock::K, base_rng::R, f::Function
  ) where {K,R}
      # Determine which RNG to use
      if crc.mode == :replaying && haskey(crc.recorder, clock)
          # ... retrieve saved RNG ...
          active_rng = saved_rng
      else
          active_rng = copy(base_rng)
          if crc.mode == :recording
              # ... save RNG snapshot ...
          end
      end

      # Call continuation with prepared RNG
      return f(active_rng)
  end

  function enable!(ctx::SamplingContext, clock, dist, te, when)
      if ctx.common_random !== nothing
          # Wrap sampler call in RNG transformation
          with_common_rng(ctx.common_random, clock, ctx.rng) do rng
              enable!(ctx.sampler, clock, dist, te, when, rng)
          end
      else
          enable!(ctx.sampler, clock, dist, te, when, ctx.rng)
      end

      # Other components use base RNG
      ctx.likelihood !== nothing && enable!(ctx.likelihood, clock, dist, te, when, ctx.rng)
      ctx.debug !== nothing && enable!(ctx.debug, clock, dist, te, when, ctx.rng)
  end
"""
mutable struct CommonRandom{K,RNG}
    record::Dict{K,Array{RNG,1}} # Value of RNG for each clock each instance.
    sample_index::Dict{K,Int} # Current number of times each clock seen.
    miss::Dict{K,Int} # Number of misses for each clock.
    mode::Symbol
    function CommonRandom{K,RNG}() where {K,RNG}
        storage = Dict{K,Array{RNG,1}}()
        index = Dict{K,Int64}()
        miss = Dict{K,Int}()
        new{K,RNG}(storage, index, miss, :record)
    end
end


clone(cr::CommonRandom{K,RNG}) where {K,RNG} = CommonRandom{K,RNG}()


"""
Doesn't reset the stored clocks, does reset miss count.
"""
function reset!(recorder::CommonRandom)
    empty!(recorder.sample_index)
    empty!(recorder.miss)
end


function reset_crn!(recorder::CommonRandom)
    reset!(recorder)
    empty!(recorder.record)
end


"""
How many times the sampler looked for a random number and found no previous value.
"""
misscount(recorder::CommonRandom) = sum(values(recorder.miss))
"""
The Pairs of clock keys and values that weren't found.
"""
misses(recorder::CommonRandom) = pairs(recorder.miss)
"""
Call this before replaying common random numbers.
Call it before each replay.
"""
function freeze_crn!(recorder::CommonRandom)
    reset!(recorder)
    recorder.mode = :replay
    return nothing
end


function with_common_rng(
    f::Function, cr::CommonRandom{K,R}, clock::K, rng::R
) where {K,R}
    clock_seen_cnt = 1 + get(cr.sample_index, clock, 0)
    samples = get(cr.record, clock, nothing)
    if samples !== nothing && clock_seen_cnt <= length(samples)
        saved_rng = copy(samples[clock_seen_cnt])
        f(saved_rng)
        # Only increment the seen count if the rng was used because Next Reaction avoids RNG.
        if saved_rng != samples[clock_seen_cnt]
            cr.sample_index[clock] = clock_seen_cnt
        end  # else don't increment.
    elseif cr.mode == :record
        rng_save = copy(rng)
        f(rng)
        if rng != rng_save
            if samples !== nothing
                push!(samples, rng_save)
            else
                cr.record[clock] = [rng_save]
            end
            cr.sample_index[clock] = clock_seen_cnt
        end
        cr.miss[clock] = get(cr.miss, clock, 0) + 1
    else  # cr.mode == :replay
        f(rng)
        cr.miss[clock] = get(cr.miss, clock, 0) + 1
    end
end
