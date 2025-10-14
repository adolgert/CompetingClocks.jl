using Base

mutable struct TrajectoryWatcher{K,T} <: EnabledWatcher{K,T}
    enabled::Dict{K,EnablingEntry{K,T}}
    loglikelihood::Float64
    curtime::Float64
    TrajectoryWatcher{K,T}() where {K,T} = new(Dict{K,EnablingEntry{K,T}}(), zero(Float64), zero(Float64))
end
export TrajectoryWatcher


function trajectoryloglikelihood(tw::TrajectoryWatcher)
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    when = tw.curtime
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


reset!(tw::TrajectoryWatcher) = (empty!(ts.enabled); tw.loglikelihood = zero(Float64); nothing)
function Base.copy!(dst::TrajectoryWatcher{K,T}, src::TrajectoryWatcher{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    dst.loglikelihood = src.loglikelihood
    return dst
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
        if haskey(ts.enabled, clock)
            delete!(ts.enabled, clock)
        end
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
        if haskey(ts.enabled, clock)
            delete!(ts.enabled, clock)
        end
    end
    ts.curtime = when
end
