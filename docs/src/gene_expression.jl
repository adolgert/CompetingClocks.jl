# # Stochastic Gene Expression
#
# This is a long example to demonstrate the use of path likelihoods to do
# importance sampling.
# The goal of the example is to estimate the probability of a transcriptional burst.
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
# ## Importance sampling
#
# The main challenge of simulating transcriptional bursts is that large bursts are rare.
# This is a good case for using importance sampling. Importance sampling is a way to simulate
# from a biased set of distributions, in this case distributions weighted towards making
# more large bursts, but to account for that bias in your final estimate.
#
# The usual way to use simulations to estimate a quantity is to sum the outcome across trajectories.
#
# ```math
# E_p[f(X)] \approx \frac{1}{N}\sum_{i=1}^N f(x_i)
# ```
# That sum is an approximation of a statistical expectation of $f$ over the probability distribution
# $p$.
# ```math
# E_p[f(X)] = \int f(x)p(x)dx
# ```
# The events in a simulation and their distributions determine the probability $p(x)$ for
# each trajectory $x$, and for rare events that probability is very small, so we will bias it.
#
# We bias it by picking different distributions for our events. Let's call the biased space
# $q$. If we work backwards from the statistical expression, we will see that we can undo
# the bias when we do our sum over trajectories.
# ```math
# E_p[f(x)] = \int f(x)\frac{p(x)}{q(x)}dx = E_q[f(X)\frac{p(X)}{q(X)}].
# ```
# Here $X$ is a sample under the distribution $q$. That means we draw from $q$ and at the
# end sum with a weight on each trajectory.
#
# ```math
# \hat{\mu} = \frac{1}{N}\sum_{i=1}^N f(X_i)w(X_i)
# ```
# Here the importance weight is $w(x)=p(x)/q(x)$. This package's samplers will calculate that
# weigth for you as a path log-likelihood, which is $\log w(x)$.
#
# Let's give this a try for gene expression. We can use an intuitive bias shift.
#
# We're going to rely on standard Julia distributions with one extra one,
# imported from the file `genedist.jl` which expresses the initial slowness
# of MRNA production using a time-varying rate.
using CompetingClocks
using Distributions
using Logging
using Random
using Printf
include("genedist.jl")
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
# You will see in the simulation that we initialize two distributions for each
# event. The first is used to determine the likelihood of the event. The second
# distribution is used to sample for the next time. These two can be the same,
# or they can differ if we want to use importance sampling.
#
# When we turn the promoter on, it turns on slowly and ramps up. This is reflected
# in a custom distribution called `TranscriptionRate`.
function turn_promoter_on(model, sampler, individual, when, θ)
    enable!(sampler, (:off, 0), [Exponential(inv(θ[:promoter_off][1])), Exponential(inv(θ[:promoter_off][2]))])
    rate1 = TranscriptionRate(θ[:transcribe_max][1], θ[:transcribe_remodel][1])
    rate2 = TranscriptionRate(θ[:transcribe_max][2], θ[:transcribe_remodel][2])
    enable!(sampler, (:transcribe, 0), [rate1, rate2])
end
#
# Turning the promoter off doesn't stop the system. We count every protein created
# until all MRNA degrade.
#
turn_promoter_off(model, sampler, individual, when, θ) = disable!(sampler, (:transcribe, 0))
#
# The transcription rate depends on the number of MRNA present as a total,
# but each MRNA can degrade individually with a Gamma distribution that starts
# it's clock when that particular MRNA is produced. The Gamma distribution is the
# same distribution we use to model how a person recovers from the flu.
#
function transcribe_mrna(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.mrna)
    total = pre_event_total
    mrna_id = length(model.mrna) + 1
    push!(model.mrna, (when, true))
    time_offset = when - 0  # The 0 is when the promoter turned on.
    rate1 = TranscriptionRate(θ[:transcribe_max][1], θ[:transcribe_remodel][1]; t0=time_offset)
    rate2 = TranscriptionRate(θ[:transcribe_max][2], θ[:transcribe_remodel][2]; t0=time_offset)
    enable!(sampler, (:transcribe, 0), [rate1, rate2])
    total = count(x -> x[2], model.mrna)
    transrate1 = Exponential(inv(θ[:translate][1] * total))
    transrate2 = Exponential(inv(θ[:translate][2] * total))
    pre_event_total > 0 && disable!(sampler, (:translate, 0))
    enable!(sampler, (:translate, 0), [transrate1, transrate2])
    gamma1 = Gamma(θ[:degrade_k][1], inv(θ[:degrade_theta][1]))
    gamma2 = Gamma(θ[:degrade_k][2], inv(θ[:degrade_theta][2]))
    enable!(sampler, (:degrade, mrna_id), [gamma1, gamma2])
end
#
# Degrading an MRNA will turn off protein production for that MRNA,
# but if it turns off the _last_ MRNA, then it turns of all translation.
#
function degrade_mrna(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.mrna)
    total = pre_event_total
    model.mrna[individual] = (zero(Time), false)
    pre_event_total > 0 && disable!(sampler, (:translate, 0))
    total = count(x -> x[2], model.mrna)
    if total > 0
        transrate1 = Exponential(inv(θ[:translate][1] * total))
        transrate2 = Exponential(inv(θ[:translate][2] * total))
        enable!(sampler, (:translate, 0), [transrate1, transrate2])
    end
end
#
# When translation happens, it isn't a one-time event. The same event
# can happen again, so we re-enable it.
#
function translate_protein(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.mrna)
    model.protein += 1
    if pre_event_total > 0
        transrate1 = Exponential(inv(θ[:translate][1] * pre_event_total))
        transrate2 = Exponential(inv(θ[:translate][2] * pre_event_total))
        enable!(sampler, (:translate, 0), [transrate1, transrate2])
    end
end
#
# Our central dynamics calls the appropriate function above when its event
# fires.
#
function step_gene!(model, sampler, which, when)
    event, individual = which
    Dict(
        :on => turn_promoter_on,
        :off => turn_promoter_off,
        :transcribe => transcribe_mrna,
        :degrade => degrade_mrna,
        :translate => translate_protein,
    )[event](model, sampler, individual, when, model.θ)
end
#
# The CompetingClocks.jl package does the sampling, but you create the mainloop.
#
function one_epoch(model, sampler)
    step_gene!(model, sampler, (:on, 0), time(sampler))
    when, which = next(sampler)
    while !isnothing(which)
        fire!(sampler, which, when)
        step_gene!(model, sampler, which, when)
        when, which = next(sampler)
    end
    basal, weighted = pathloglikelihood(sampler, time(sampler))
    logimportance = basal - weighted
    return (model.protein, logimportance)
end
#
# This runs the simulation many times with chosen parameters.
#
function run_epochs(epoch_cnt, use_importance, rng)
    ## We define two sets of parameters. The first biases the simulation towards
    ## producing a rare event and the second is the basal rate we use to evaluate
    ## the importance of those events.
    params = Dict(
        :promoter_off => (0.6, 0.2), # per minute
        :transcribe_max => (10.0, 10.0), # mRNA/min
        :transcribe_remodel => (1.0, 1.0), # per minute, rate of chromatin opening.
        :degrade_k => (4.0, 4.0),  # k for Gamma
        :degrade_theta => (4 * 4 / 30, 4 * 4 / 30), # theta for Gamma
        :translate => (1.00, 1.05), # proteins/min/mRNA
    )
    protein = zeros(Int, epoch_cnt)
    importance = zeros(Float64, epoch_cnt)
    model = GeneExpression(params)
    builder = SamplerBuilder(
        Tuple{Symbol,Int}, Float64;
        sampler_spec=:firsttofire,
        path_likelihood=true,
        likelihood_cnt=2,
    )
    sampler = SamplingContext(builder, rng)
    sample_from_distribution!(sampler, use_importance ? 2 : 1)
    for epoch_idx in eachindex(protein)
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
#
# Let's do one run without importance sampling and print the results to see how
#  often the rare event happens. We'll see it's fairly rare.
#
with_logger(ConsoleLogger(stdout, Logging.Info)) do
    observed, importance = run_epochs(100, false, Xoshiro(324923))
    show_observed(observed)
end
#
# If we run it ten times with importance sampling, we can see how good the
# statistics get. I'm keepign the sample count low here because this is just
# documentation.
#
function variations(var_cnt, N)
    prob_over_1000 = zeros(Float64, var_cnt)
    fraction_over = zeros(Float64, var_cnt)
    rng = Xoshiro(234291022)
    for pidx in eachindex(prob_over_1000)
        observed, Δ = run_epochs(N, true, rng)
        ## Use log-space trick to avoid summing a bunch of zeros and extremely small numbers.
        importance = exp.(Δ .- maximum(Δ))
        ## This is the self-normalized estimator.
        prob_over_1000[pidx] = sum((observed .>= 1000) .* importance) / sum(importance)
        ## The unbiased estimator uses 1/N.
        ## prob_over_1000[pidx] = sum((observed .>= 1000) .* importance) / N
        fraction_over[pidx] = count(x -> x > 1000, observed) / length(observed)
        println("mean weight $(mean(importance))")
        println("ESS $(sum(importance)^2 / sum(importance.^2))")
    end
    println("fraction_over")
    println(join([@sprintf("%.2g", x) for x in fraction_over], ", "))
    println("probability_over")
    println(join([@sprintf("%.2g", x) for x in prob_over_1000], ", "))
end
variations(10, 1_000)
#
# ## References
#
# * Zong, Chenghang, Lok‐hang So, Leonardo A Sepúlveda, Samuel O Skinner, and Ido Golding. “Lysogen Stability Is Determined by the Frequency of Activity Bursts from the Fate‐determining Gene.” Molecular Systems Biology 6, no. 1 (2010): 440. https://doi.org/10.1038/msb.2010.96.
# * Raj, Arjun, and Alexander Van Oudenaarden. “Nature, Nurture, or Chance: Stochastic Gene Expression and Its Consequences.” Cell 135, no. 2 (2008): 216–26. https://doi.org/10.1016/j.cell.2008.09.050.
# * Cai, Long, Nir Friedman, and X. Sunney Xie. “Stochastic Protein Expression in Individual Cells at the Single Molecule Level.” Nature 440, no. 7082 (2006): 358–62. https://doi.org/10.1038/nature04599.
# * Horowitz, Jordan M, and Rahul V Kulkarni. “Stochastic Gene Expression Conditioned on Large Deviations.” Physical Biology 14, no. 3 (2017): 03LT01. https://doi.org/10.1088/1478-3975/aa6d89.
# * McAdams, Harley H., and Adam Arkin. “Stochastic Mechanisms in Gene Expression.” Proceedings of the National Academy of Sciences 94, no. 3 (1997): 814–19. https://doi.org/10.1073/pnas.94.3.814.
