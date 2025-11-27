using StatsBase


struct SamplerSUT
    travel::TravelConfig
    state_cnt::Int
    step_cnt::Int
end


mutable struct ScoreEvents{T}
    tracker::T
    last_state::Int
    last_time::Float64
    # This should be uniformly distributed.
    waiting_cumulant::Dict{Int,Vector{Float64}}
    step_cnt::Int
    mark_brier::Float64
    mark_null::Vector{Float64}
    rng::Xoshiro
    ScoreEvents(tracker::TR, rng) where {TR} = new{TR}(
        tracker, 0, 0.0, Dict{Int,Vector{Float64}}(), 0, 0.0, zeros(1999), rng
        )
end


function (score::ScoreEvents)(when::Float64, which::Int64)
    # The first part is the Doob-Meyer calculation for waiting times.
    score.step_cnt += 1
    transition = score.last_state
    cumulant = CompetingClocks.stepcumulant(score.tracker, score.last_time, when)
    if haskey(score.waiting_cumulant, transition)
        push!(score.waiting_cumulant[transition], cumulant)
    else
        score.waiting_cumulant[transition] = [cumulant]
    end

    probs = CompetingClocks.stepconditionalprobability(score.tracker, when)
    # Sample these 1999 times. Calculate Brier contributions.
    briers = Dict(
        k => sum(cp -> (cp[2] - (cp[1] == k ? 1 : 0))^2, probs)
        for k in keys(probs)
    )
    ks = collect(keys(probs))
    score.mark_brier += briers[which]
    score.mark_null .+= sample(score.rng, [briers[k] for k in ks], Weights([probs[j] for j in ks]))

    score.last_state = which
    score.last_time = when
end


"""
This is a Doob-Meyer test on the waiting time out of each state.
"""
function waiting_metric(score::ScoreEvents)
    measures = Any[]
    for idx in keys(score.waiting_cumulant)
        test = ApproximateOneSampleKSTest(score.waiting_cumulant[idx], Uniform(0, 1))
        p = pvalue(test)
        push!(measures, (; pvalue=p, clock=idx, count = length(score.waiting_cumulant[idx]), supremum_of_difference=test.δ, test=test))
    end
    return measures
end


function mark_calibration(score::ScoreEvents)
    B = length(score.mark_null)
    return (;
        pvalue=(1 + count(x -> x >= score.mark_brier, score.mark_null)) / (B + 1),
        count = score.step_cnt,
        clock = 0
    )
end


function collect_data_single(smethod, sut::SamplerSUT, rng)
    sampler = smethod(Int,Float64)
    model = Travel(sut.state_cnt, sut.travel, rng)
    tracker = TrackWatcher{Int64,Float64}()
    observer = ScoreEvents(tracker, rng)
    recorder = SamplerRecord(tracker, rng)
    travel_run(10_000, sampler, model, observer, recorder, rng)
    metrics = waiting_metric(observer)
    push!(metrics, mark_calibration(observer))
    return metrics
end
