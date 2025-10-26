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
# ## State of the system
#
# We start the system when the promoter turns ON and stop the simulation when
# the promoter turns OFF.
#
#  1. Vector of MRNA. Each one was created at a certain time and is enabled/disabled.
#  2. Count of total proteins created.
#
# ## Events in the system
#
#  1. `(:on, 0)`, Turn on the promoter. We use this to start the simulation.
#  2. `(:transcribe, 0)` - When this fires, the promoter creates an MRNA.
#  3. `(:translate, 0)` - The rate of translation is proportional to the number
#     of MRNA that currently exist.
#  4. `(:degrade, mrna_id)` - A particular MRNA will degrade, turning it off.
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
    θ::Dict{Symbol,NTuple{2,Float64}}
    function GeneExpression(params)
        mrna = Tuple{Time,Bool}[]
        sizehint!(mrna, 2000)
        new(mrna, 0, params)
    end
end
Base.empty!(ge::GeneExpression) = (empty!(ge.mrna); ge.protein = 0; nothing)
#
# You will see in the simulation that we initialize two distributions. The first
# distribution is used to sample for the next time. The second distribution
# is used to determine the likelihood of the event. These two can be the same,
# or they can differ if we want to use importance sampling.
#
function step_gene!(model, sampler, which, when)
    θ = model.θ
    event, individual = which
    if event == :on
        enable!(sampler, (:off, 0), [Exponential(inv(θ[:promoter_off][1])), Exponential(inv(θ[:promoter_off][2]))])
        rate1 = TranscriptionRate(θ[:transcribe_max][1], θ[:transcribe_remodel][1])
        rate2 = TranscriptionRate(θ[:transcribe_max][2], θ[:transcribe_remodel][2])
        enable!(sampler, (:transcribe, 0), [rate1, rate2])
    elseif event == :off
        disable!(sampler, (:transcribe, 0))
    elseif event == :transcribe
        mrna_id = length(model.mrna) + 1
        push!(model.mrna, (when, true))
        time_offset = when - 0  # The 0 is when the promoter turned on.
        rate1 = TranscriptionRate(θ[:transcribe_max][1], θ[:transcribe_remodel][1]; t0=time_offset)
        rate2 = TranscriptionRate(θ[:transcribe_max][2], θ[:transcribe_remodel][2]; t0=time_offset)
        enable!(sampler, (:transcribe, 0), [rate1, rate2])
        total = count(x -> x[2], model.mrna)
        transrate1 = Exponential(inv(θ[:translate][1] * total))
        transrate2 = Exponential(inv(θ[:translate][2] * total))
        enable!(sampler, (:translate, 0), [transrate1, transrate2])
        # Julia uses shape and scale, but we specify rate, so use 1-over.
        gamma1 = Gamma(θ[:degrade_k][1], inv(θ[:degrade_theta][1]))
        gamma2 = Gamma(θ[:degrade_k][2], inv(θ[:degrade_theta][2]))
        enable!(sampler, (:degrade, mrna_id), [gamma1, gamma2])
    elseif event == :degrade
        model.mrna[individual] = (zero(Time), false)
        total = count(x -> x[2], model.mrna)
        if total > 0
            transrate1 = Exponential(inv(θ[:translate][1] * total))
            transrate2 = Exponential(inv(θ[:translate][2] * total))
            enable!(sampler, (:translate, 0), [transrate1, transrate2])
        end
    elseif event == :translate
        model.protein += 1
    end
end

function one_epoch(model, sampler)
    step_gene!(model, sampler, (:on, 0), time(sampler))
    when, which = next(sampler)
    while !isnothing(which)
        fire!(sampler, which, when)
        step_gene!(model, sampler, which, when)
        when, which = next(sampler)
    end
    # The first is the one we sampled. The second is the basal rates.
    weighted, basal = trajectoryloglikelihood(sampler, when)
    importance = exp(basal - weighted)
    return (model.protein, importance)
end

function run_epochs(epoch_cnt, importance, rng)
    # We define two sets of parameters. The first biases the simulation towards
    # producing a rare event and the second is the basal rate we use to evaluate
    # the importance of those events.
    params = Dict(
        :promoter_off => (0.1, 0.2), # per minute
        :transcribe_max => (10.0, 10.0), # mRNA/min
        :transcribe_remodel => (1.0, 1.0), # per minute, rate of chromatin opening.
        :degrade_k => (4.0, 4.0),  # k for Gamma
        :degrade_theta => (4 / 30, 4 / 30), # theta for Gamma
        :translate => (2.0, 2.0), # proteins/min/mRNA
    )
    if !importance
        # Erase the weighted params
        params = Dict((k => (v[2], v[2])) for (k, v) in params)
        println("erasing importance")
    end
    protein = zeros(Int, epoch_cnt)
    importance = zeros(Float64, epoch_cnt)
    for epoch_idx in eachindex(protein)
        model = GeneExpression(params)
        builder = SamplerBuilder(
            Tuple{Symbol,Int}, Float64;
            sampler_spec=:firsttofire,
            trajectory_likelihood=true,
            likelihood_cnt=2,
        )
        sampler = SamplingContext(builder, rng)
        (cnt, weight) = one_epoch(model, sampler)
        protein[epoch_idx] = cnt
        importance[epoch_idx] = weight
        empty!(model)
        reset!(sampler)
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
observed, importance = run_epochs(1_000, false, Xoshiro(324923))
show_observed(observed)

function variations(var_cnt)
    prob_over_1000 = zeros(Float64, var_cnt)
    fraction_over = zeros(Float64, var_cnt)
    rng = Xoshiro(234291022)
    for pidx in eachindex(prob_over_1000)
        observed, importance = run_epochs(10000, true, rng)
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
