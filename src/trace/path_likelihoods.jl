using Base

struct PathEntry{K,T}
    clock::K
    distribution::Vector{UnivariateDistribution}
    te::T    # The zero-point, in absolute time, for the distribution.
    when::T  # When the distribution was enabled.
end


mutable struct PathLikelihoods{K,T}
    enabled::Dict{K,PathEntry{K,T}}
    loglikelihood::Vector{Float64}
    curtime::Float64
    PathLikelihoods{K,T}(cnt) where {K,T} = new(Dict{K,PathEntry{K,T}}(), zeros(Float64, cnt), zero(Float64))
end
export PathLikelihoods

_likelihood_cnt(pl::PathLikelihoods) = length(pl.loglikelihood)

Base.length(pl::PathLikelihoods) = length(pl.enabled)
enabled(pl::PathLikelihoods) = keys(pl.enabled)
isenabled(pl::PathLikelihoods, clock) = haskey(pl.enabled, clock)

function Base.copy!(dst::PathLikelihoods{K,T}, src::PathLikelihoods{K,T}) where {K,T}
    copy!(dst.enabled, src.enabled)
    copy!(dst.loglikelihood, src.loglikelihood)
    dst.curtime = src.curtime
    return dst
end


function trajectoryloglikelihood(tw::PathLikelihoods, when)
    # When this is called, there will be transitions that have not yet fired, and
    # they need to be included as though they were just disabled.
    if when > tw.curtime
        remaining = zeros(Float64, _likelihood_cnt(tw))
        for entry in values(tw.enabled)
            for idx in eachindex(remaining)
                if when > entry.te
                    remaining[idx] += logccdf(entry.distribution[idx], when - entry.te)
                    if entry.when > entry.te
                        remaining[idx] -= logccdf(entry.distribution[idx], entry.when - entry.te)
                    end
                end
            end
        end
        return tw.loglikelihood .+ remaining
    else
        return tw.loglikelihood
    end
end


function reset!(tw::PathLikelihoods)
    empty!(tw.enabled)
    tw.loglikelihood .= zero(Float64)
    tw.curtime = 0.0
    nothing
end


function enable!(ts::PathLikelihoods{K,T}, clock::K, dist::UnivariateDistribution, te::T, when::T, rng::AbstractRNG) where {K,T}
    ts.enabled[clock] = PathEntry{K,T}(clock, UnivariateDistribution[copy(dist)], te, when)
end


function enable!(ts::PathLikelihoods{K,T}, clock::K, dist::Vector, te::T, when::T, rng::AbstractRNG) where {K,T}
    ts.enabled[clock] = PathEntry{K,T}(clock, copy(dist), te, when)
end


function disable!(ts::PathLikelihoods{K,T}, clock::K, when::T) where {K,T}
    entry = get(ts.enabled, clock, nothing)
    isnothing(entry) && return
    if length(entry.distribution) == 1
        log_delta = zero(Float64)
        if when > entry.te
            log_delta += logccdf(entry.distribution[1], when - entry.te)
            if entry.when > entry.te
                log_delta -= logccdf(entry.distribution[1], entry.when - entry.te)
            end
        end
        ts.loglikelihood .+= log_delta
    else
        for idx in eachindex(entry.distribution)
            if when > entry.te
                ts.loglikelihood[idx] += logccdf(entry.distribution[idx], when - entry.te)
                if entry.when > entry.te
                    ts.loglikelihood[idx] -= logccdf(entry.distribution[idx], entry.when - entry.te)
                end
            end
        end
    end
    delete!(ts.enabled, clock)
end


function fire!(ts::PathLikelihoods{K,T}, clock::K, when::T) where {K,T}
    entry = get(ts.enabled, clock, nothing)
    if !isnothing(entry)
        if length(entry.distribution) == 1
            log_delta = zero(Float64)
            if when > entry.te
                log_delta += logpdf(entry.distribution[1], when - entry.te)
            end
            # Adjust for an enabling time that was shifted left.
            if entry.when > entry.te
                log_delta -= logccdf(entry.distribution[1], entry.when - entry.te)
            end
            ts.loglikelihood .+= log_delta
        else
            for idx in eachindex(entry.distribution)
                if when > entry.te
                    ts.loglikelihood[idx] += logpdf(entry.distribution[idx], when - entry.te)
                end
                # Adjust for an enabling time that was shifted left.
                if entry.when > entry.te
                    ts.loglikelihood[idx] -= logccdf(entry.distribution[idx], entry.when - entry.te)
                end
            end
        end
        delete!(ts.enabled, clock)
    end
    ts.curtime = when
end
