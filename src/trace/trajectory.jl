using Base

"""
    TrajectoryWatcher{K,T,L}

This doesn't sample but calculates the likelihood of the path of samples from
start to finish. It has many of the same interface functions as a sampler,
but the core value is in the [`pathloglikelihood`](@ref) function.

The `L` type parameter is the number type of the accumulated log-likelihood. It
defaults to `Float64`, but allowing it to vary is what lets `ForwardDiff.Dual`
values from differentiated distribution parameters flow through the accumulator.
"""
mutable struct TrajectoryWatcher{K,T,L} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    loglikelihood::L
    curtime::T
    TrajectoryWatcher{K,T,L}() where {K,T,L} = new(Dict{K,EnablingEntry{K,T}}(), zero(L), zero(T))
end

# Back-compat: existing call sites construct without specifying the accumulator
# type, which should keep meaning a Float64 accumulator.
TrajectoryWatcher{K,T}() where {K,T} = TrajectoryWatcher{K,T,Float64}()


clone(tw::TrajectoryWatcher{K,T,L}) where {K,T,L} = TrajectoryWatcher{K,T,L}()


function copy_clocks!(dst::TrajectoryWatcher{K,T,L}, src::TrajectoryWatcher{K,T,L}) where {K,T,L}
    copy!(dst.enabled, src.enabled)
    dst.loglikelihood = src.loglikelihood
    dst.curtime = src.curtime
    return dst
end


function pathloglikelihood(tw::TrajectoryWatcher{K,T,L}, when) where {K,T,L}
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    remaining = zero(L)
    for entry in values(tw.enabled)
        if when > entry.te
            remaining += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                remaining -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
    end

    return tw.loglikelihood + remaining
end


function reset!(tw::TrajectoryWatcher{K,T,L}) where {K,T,L}
    empty!(tw.enabled)
    tw.loglikelihood = zero(L)
    tw.curtime = zero(T)
    nothing
end

function disable!(ts::TrajectoryWatcher{K,T}, clock::K, when::T) where {K,T}
    entry = get(ts.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
        delete!(ts.enabled, clock)
    else
        error("Cannot disable clock $clock at $when because it is not enabled.")
    end
end


function fire!(ts::TrajectoryWatcher{K,T}, clock::K, when::T) where {K,T}
    entry = get(ts.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logpdf(entry.distribution, when - entry.te)
        end
        # Adjust for an enabling time that was shifted left.
        if entry.when > entry.te
            ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
        end
        delete!(ts.enabled, clock)
    else
        error("Cannot fire clock $clock at $when because it is not enabled.")
    end
    ts.curtime = when
end
