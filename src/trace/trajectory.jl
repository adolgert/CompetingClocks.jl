using Base
using Distributions: UnivariateDistribution

mutable struct TrajectoryWatcher{K,T}
    track::TrackWatcher{K,T}
    loglikelihood::Float64
    curtime::Float64
    TrajectoryWatcher{K,T}() where {K,T} = new(TrackWatcher{K,T}(), zero(Float64), zero(Float64))
end
export TrajectoryWatcher


function trajectoryloglikelihood(tw::TrajectoryWatcher)
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    when = tw.curtime
    remaining = zero(Float64)
    for entry in values(tw.track.enabled)
        if when > entry.te
            remaining += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                remaining -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
    end

    return tw.loglikelihood + remaining
end


reset!(tw::TrajectoryWatcher) = (reset!(tw.track); tw.loglikelihood=zero(Float64); nothing)
function Base.copy!(dst::TrajectoryWatcher{K,T}, src::TrajectoryWatcher{K,T}) where {K,T}
    copy!(dst.track, src.track)
    dst.loglikelihood = src.loglikelihood
end

Base.iterate(ts::TrajectoryWatcher) = iterate(ts.track)
Base.iterate(ts::TrajectoryWatcher, i::Int64) = iterate(ts.track, i)
Base.length(ts::TrajectoryWatcher) = length(ts.track)


function enable!(ts::TrajectoryWatcher{K,T}, clock::K, dist::UnivariateDistribution, te::T, when::T, rng) where {K,T}
    enable!(ts.track, clock, dist, te, when, rng)
end


function disable!(ts::TrajectoryWatcher{K,T}, clock::K, when) where {K,T}
    entry = get(ts.track.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logccdf(entry.distribution, when - entry.te)
            if entry.when > entry.te
                ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
            end
        end
        disable!(ts.track, clock, when)
    end
end


function fire!(ts::TrajectoryWatcher{K,T}, clock::K, when) where {K,T}
    entry = get(ts.track.enabled, clock, nothing)
    if entry !== nothing
        if when > entry.te
            ts.loglikelihood += logpdf(entry.distribution, when - entry.te)
        end
        # Adjust for an enabling time that was shifted left.
        if entry.when > entry.te
            ts.loglikelihood -= logccdf(entry.distribution, entry.when - entry.te)
        end
        disable!(ts.track, clock, when)
    end
    ts.curtime = when
end


function steploglikelihood(tw::TrajectoryWatcher{K,T}, now, when, which_fires) where {K,T}
    steploglikelihood(tw.track, now, when, which_fires)
end
