# # Stochastic Gene Expression
#
# The goal is to estimate the probability of a transcriptional burst.
#
# DNA becomes proteins in a three-step process.
#
#  1. DNA can replicate itself.
#  2. Transcription changes a gene in DNA into RNA.
#  3. Translation turns RNA into protein.
#
# We will account for the changing states of the system.
#
#  * A **gene** can be ON or OFF. When it's ON, it can make RNA.
#    We are looking for the case where the switch stays on longer.
#
#  * The gene has startup time to transcribe its first RNA. One possibile
#    hazard rate comes from the LogNormal distribution.
#
#  * mRNA degrades. Different parts of the mRNA are degrading individually with
#    constant rates until they reach a threshold. The time to the threshold
#    is described by a Gamma distribution.
#
#   * Each mRNA makes proteins until it degrades.
#
# State of the system:
#
using CompetingClocks
using Distributions
using Random
Time = Float64
Epoch = Int
mutable struct GeneExpression
    promoter::Bool
    promoter_on_time::Time
    promoter_epoch::Epoch
    mrna::Dict{Epoch,Vector{Tuple{Time,Bool}}}
    θ::Dict{Symbol,Float64}
    GeneExpression(params) = new(false, zero(Time), 0, Dict{Epoch,Vector{Tuple{Time,Bool}}}(), params)
end
mutable struct GeneObserver
    protein::Vector{Int}
    completed::Epoch
    GeneObserver() = new(Int[], 0)
end
#
function step_gene!(model, sampler, which, when)
    epoch = model.promoter_epoch
    θ = model.θ
    event, event_epoch, individual = which
    if event == :on
        model.promoter = true
        model.promoter_epoch += 1
        model.promoter_on_time = when
        enable!(sampler, (:off, model.promoter_epoch, 0), Exponential(inv(θ[:promoter_off])))
        rate = TranscriptionRate(θ[:transcribe_max], θ[:transcribe_remodel])
        enable!(sampler, (:transcribe, model.promoter_epoch, 0), rate)
    elseif event == :off
        model.promoter = false
        enable!(sampler, (:on, model.promoter_epoch, 0), Exponential(inv(θ[:promoter_on])))
        disable!(sampler, (:transcribe, epoch, 0))
    elseif event == :transcribe
        times = get!(model.mrna, epoch, Tuple{Time,Bool}[])
        mrna_id = length(times) + 1
        push!(times, (when, true))
        time_offset = when - model.promoter_on_time
        rate = TranscriptionRate(θ[:transcribe_max], θ[:transcribe_remodel]; t0=time_offset)
        enable!(sampler, (:transcribe, model.promoter_epoch, 0), rate)
        total = count(x -> x[2], model.mrna[epoch])
        transrate = Exponential(inv(θ[:translate] * total))
        enable!(sampler, (:translate, epoch, 0), transrate)
        # Julia uses shape and scale, but we specify rate, so use 1-over.
        gamma = Gamma(θ[:degrade_k], inv(θ[:degrade_theta]))
        enable!(sampler, (:degrade, epoch, mrna_id), gamma)
    elseif event == :degrade
        model.mrna[event_epoch][individual] = (zero(Time), false)
        total = count(x -> x[2], model.mrna[event_epoch])
        if total > 0
            transrate = Exponential(inv(θ[:translate] * total))
            enable!(sampler, (:translate, event_epoch, 0), transrate)
        else
            delete!(model.mrna, event_epoch)
        end
    elseif event == :translate
        # Here the second argument to the `which` is the epoch
    end
end
function observe_gene!(model, observer, which, when)
    event, event_epoch, individual = which
    if event == :on
        push!(observer.protein, 0)
    elseif event == :degrade && !haskey(model.mrna, event_epoch)
        observer.completed += 1
        return observer.completed
    elseif event == :translate
        observer.protein[event_epoch] += 1
    end
    return nothing
end
function run_epochs(epoch_cnt)
    params = Dict(
        :promoter_on => 0.04, # per minute
        :promoter_off => 0.2, # per minute
        :promoter_off_biased => 0.02, # per minute
        :transcribe_max => 10.0, # mRNA/min
        :transcribe_remodel => 1.0, # per minute, rate of chromatin opening.
        :degrade_k => 4,  # k for Gamma
        :degrade_theta => 4 / 30, # theta for Gamma
        :translate => 2, # proteins/min/mRNA
    )
    model = GeneExpression(params)
    observer = GeneObserver()
    builder = SamplerBuilder(Tuple{Symbol,Int,Int}, Float64; sampler_spec=:firsttofire)
    rng = Xoshiro(92734924)
    sampler = SamplingContext(builder, rng)
    step_gene!(model, sampler, (:on, 0, 0), time(sampler))
    observe_gene!(model, observer, (:on, 0, 0), 0.0)
    running = true
    while running
        when, which = next(sampler)
        fire!(sampler, which, when)
        step_gene!(model, sampler, which, when)
        running = observe_gene!(model, observer, which, when) != epoch_cnt
    end
    return observer
end
function show_observed(observer)
    bins = 100 * collect(1:10)
    keep = observer.protein[1:observer.completed]
    gt_bin = [sum(keep .> bin) for bin in bins]
    for idx in eachindex(bins)
        println("bin $(bins[idx]) count $(gt_bin[idx])")
    end
    println("total $(observer.completed)")
end
observed = run_epochs(100 * 300)
show_observed(observed)
