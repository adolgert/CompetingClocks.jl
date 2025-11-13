using Base

"""
    TrajectoryWatcher{K,T}

This doesn't sample but calculates the likelihood of the path of samples from
start to finish. It has many of the same interface functions as a sampler,
but the core value is in the [`pathloglikelihood`](@ref) function.
"""
mutable struct TrajectoryWatcher{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    loglikelihood::Float64
    curtime::T
    TrajectoryWatcher{K,T}() where {K,T} = new(Dict{K,EnablingEntry{K,T}}(), zero(Float64), zero(T))
end
export TrajectoryWatcher


clone(tw::TrajectoryWatcher{K,T}) where {K,T} = TrajectoryWatcher{K,T}()


function copy_clocks!(dst::TrajectoryWatcher{K,T}, src::TrajectoryWatcher{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    dst.loglikelihood = src.loglikelihood
    dst.curtime = src.curtime
    return dst
end


function pathloglikelihood(tw::TrajectoryWatcher, when)
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    remaining = zero(Float64)
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


function reset!(tw::TrajectoryWatcher{K,T}) where {K,T}
    empty!(tw.enabled)
    tw.loglikelihood = zero(Float64)
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
