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
using Printf
Time = Float64
Epoch = Int
mutable struct GeneExpression
    mrna::Vector{Tuple{Time,Bool}}
    protein::Int
    θ::Dict{Symbol,Float64}
    function GeneExpression(params)
        mrna = Tuple{Time,Bool}[]
        sizehint!(mrna, 2000)
        new(mrna, 0, params)
    end
end
Base.empty!(ge::GeneExpression) = (empty!(ge.mrna); ge.protein = 0; nothing)
#
function step_gene!(model, sampler, which, when)
    θ = model.θ
    event, individual = which
    if event == :on
        enable!(sampler, (:off, 0), Exponential(inv(θ[:promoter_off])))
        rate = TranscriptionRate(θ[:transcribe_max], θ[:transcribe_remodel])
        enable!(sampler, (:transcribe, 0), rate)
    elseif event == :off
        disable!(sampler, (:transcribe, 0))
    elseif event == :transcribe
        mrna_id = length(model.mrna) + 1
        push!(model.mrna, (when, true))
        time_offset = when - 0  # The 0 is when the promoter turned on.
        rate = TranscriptionRate(θ[:transcribe_max], θ[:transcribe_remodel]; t0=time_offset)
        enable!(sampler, (:transcribe, 0), rate)
        total = count(x -> x[2], model.mrna)
        transrate = Exponential(inv(θ[:translate] * total))
        enable!(sampler, (:translate, 0), transrate)
        # Julia uses shape and scale, but we specify rate, so use 1-over.
        gamma = Gamma(θ[:degrade_k], inv(θ[:degrade_theta]))
        enable!(sampler, (:degrade, mrna_id), gamma)
    elseif event == :degrade
        model.mrna[individual] = (zero(Time), false)
        total = count(x -> x[2], model.mrna)
        if total > 0
            transrate = Exponential(inv(θ[:translate] * total))
            enable!(sampler, (:translate, 0), transrate)
        end
    elseif event == :translate
        model.protein += 1
    end
end

function one_epoch(model, sampler, model_weighted, sampler_weighted)
    step_gene!(model, sampler, (:on, 0), time(sampler))
    step_gene!(model_weighted, sampler_weighted, (:on, 0), time(sampler_weighted))
    when, which = next(sampler_weighted)
    while !isnothing(which)
        fire!(sampler_weighted, which, when)
        fire!(sampler, which, when)
        step_gene!(model_weighted, sampler_weighted, which, when)
        step_gene!(model, sampler, which, when)
        when, which = next(sampler_weighted)
    end
    weighted = trajectoryloglikelihood(sampler_weighted, when)
    basal = trajectoryloglikelihood(sampler, when)
    importance = exp(basal - weighted)
    return (model.protein, importance)
end

function run_epochs(epoch_cnt, rng)
    params = Dict(
        :promoter_off => 0.2, # per minute
        :transcribe_max => 10.0, # mRNA/min
        :transcribe_remodel => 1.0, # per minute, rate of chromatin opening.
        :degrade_k => 4,  # k for Gamma
        :degrade_theta => 4 / 30, # theta for Gamma
        :translate => 2, # proteins/min/mRNA
    )
    weighted_params = Dict(
        :promoter_off => 0.05, # per minute
        :transcribe_max => 15.0, # mRNA/min
        :transcribe_remodel => 1.0, # per minute, rate of chromatin opening.
        :degrade_k => 4,  # k for Gamma
        :degrade_theta => 4 / 30, # theta for Gamma
        :translate => 2, # proteins/min/mRNA
    )
    model = GeneExpression(params)
    model_weighted = GeneExpression(weighted_params)
    builder = SamplerBuilder(
        Tuple{Symbol,Int}, Float64;
        sampler_spec=:firsttofire,
        trajectory_likelihood=true,
    )
    sampler = SamplingContext(builder, rng)
    sampler_weighted = SamplingContext(builder, rng)
    protein = zeros(Int, epoch_cnt)
    importance = zeros(Float64, epoch_cnt)
    for epoch_idx in eachindex(protein)
        (cnt, weight) = one_epoch(model, sampler, model_weighted, sampler_weighted)
        protein[epoch_idx] = cnt
        importance[epoch_idx] = weight
        empty!(model)
        empty!(model_weighted)
        reset!(sampler)
        reset!(sampler_weighted)
    end
    return protein, importance
end
function show_observed(observed)
    bins = 100 * collect(1:10)
    gt_bin = [sum(observed .> bin) for bin in bins]
    for idx in eachindex(bins)
        println("bin $(bins[idx]) count $(gt_bin[idx])")
    end
    println("total $(length(observed))")
end
observed, importance = run_epochs(100, Xoshiro(324923))
show_observed(observed)

function variations(var_cnt)
    prob_over_1000 = zeros(Float64, var_cnt)
    fraction_over = zeros(Float64, var_cnt)
    rng = Xoshiro(234291022)
    for pidx in eachindex(prob_over_1000)
        observed, importance = run_epochs(10000, rng)
        # This is the self-normalized estimator. The unbiased estimator uses 1/N.
        prob_over_1000[pidx] = sum((observed .>= 1000) .* importance) / sum(importance)
        fraction_over[pidx] = count(x -> x > 1000, observed) / length(observed)
    end
    println("fraction_over")
    println(join([@sprintf("%.2g", x) for x in fraction_over], ", "))
    println("probability_over")
    println(join([@sprintf("%.2g", x) for x in prob_over_1000], ", "))
end
variations(10)
#
# ## References
#
# * Zong, Chenghang, Lok‐hang So, Leonardo A Sepúlveda, Samuel O Skinner, and Ido Golding. “Lysogen Stability Is Determined by the Frequency of Activity Bursts from the Fate‐determining Gene.” Molecular Systems Biology 6, no. 1 (2010): 440. https://doi.org/10.1038/msb.2010.96.
# * Raj, Arjun, and Alexander Van Oudenaarden. “Nature, Nurture, or Chance: Stochastic Gene Expression and Its Consequences.” Cell 135, no. 2 (2008): 216–26. https://doi.org/10.1016/j.cell.2008.09.050.
# * Cai, Long, Nir Friedman, and X. Sunney Xie. “Stochastic Protein Expression in Individual Cells at the Single Molecule Level.” Nature 440, no. 7082 (2006): 358–62. https://doi.org/10.1038/nature04599.
# * Horowitz, Jordan M, and Rahul V Kulkarni. “Stochastic Gene Expression Conditioned on Large Deviations.” Physical Biology 14, no. 3 (2017): 03LT01. https://doi.org/10.1088/1478-3975/aa6d89.
# * McAdams, Harley H., and Adam Arkin. “Stochastic Mechanisms in Gene Expression.” Proceedings of the National Academy of Sciences 94, no. 3 (1997): 814–19. https://doi.org/10.1073/pnas.94.3.814.
