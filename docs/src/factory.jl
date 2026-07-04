# # Assembly Line
#
# This is a long example to demonstrate the use of path likelihoods to do
# importance sampling.
# The goal of the example is to estimate the probability of a production burst.
#
# Raw parts become finished widgets in a two-step production pipeline.
#
#  1. A raw machine fabricates parts.
#  2. Parts are assembled into widgets.
#
# We will account for the changing states of the system.
#
#  * A **machine** can be RUNNING or IDLE. When it's running, it can make parts.
#    We are looking for the case where the switch stays on longer.
#
#  * The machine has a startup time to fabricate its first part. One possible
#    hazard rate comes from the LogNormal distribution.
#
#  * Parts wear out. Different pieces of a part are wearing out individually with
#    constant rates until they reach a threshold. The time to the threshold
#    is described by a Gamma distribution.
#
#   * Each part is assembled into widgets until it is scrapped.
#
# ## Importance sampling
#
# The main challenge of simulating production bursts is that large bursts are rare.
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
# Let's give this a try for a workstation. We can use an intuitive bias shift.
#
# We're going to rely on standard Julia distributions with one extra one,
# imported from the file `setupdist.jl` which expresses the initial slowness
# of part production using a time-varying rate.
using CompetingClocks
using Distributions
using Logging
using Random
using Printf
include("setupdist.jl")
#
# ## State of the system
#
# We start the system when the machine turns ON and stop the simulation when
# the machine turns OFF.
#
#  1. Vector of parts. Each one was created at a certain time and is enabled/disabled.
#  2. Count of total widgets created.
#
# ## Events in the system
#
#  1. `(:on, 0)`, Turn on the machine. We use this to start the simulation.
#  2. `(:fabricate, 0)` - When this fires, the machine creates a part.
#  3. `(:assemble, 0)` - The rate of assembly is proportional to the number
#     of parts that currently exist.
#  4. `(:scrap, part_id)` - A particular part will wear out, turning it off.
#
Time = Float64
Epoch = Int
mutable struct Workstation
    parts::Vector{Tuple{Time,Bool}}
    widgets::Int
    θ::Dict{Symbol,NTuple{2,Float64}}
    function Workstation(params)
        parts = Tuple{Time,Bool}[]
        sizehint!(parts, 2000)
        new(parts, 0, params)
    end
end
Base.empty!(ws::Workstation) = (empty!(ws.parts); ws.widgets = 0; nothing)
#
# You will see in the simulation that we initialize two distributions for each
# event. The first is used to determine the likelihood of the event. The second
# distribution is used to sample for the next time. These two can be the same,
# or they can differ if we want to use importance sampling.
#
# When we start the machine, it turns on slowly and ramps up. This is reflected
# in a custom distribution called `SetupRate`.
function start_machine(model, sampler, individual, when, θ)
    enable!(sampler, (:off, 0), [Exponential(inv(θ[:machine_off][1])), Exponential(inv(θ[:machine_off][2]))])
    rate1 = SetupRate(θ[:fabricate_max][1], θ[:fabricate_setup][1])
    rate2 = SetupRate(θ[:fabricate_max][2], θ[:fabricate_setup][2])
    enable!(sampler, (:fabricate, 0), [rate1, rate2])
end
#
# Turning the machine off doesn't stop the system. We count every widget created
# until all parts are scrapped.
#
stop_machine(model, sampler, individual, when, θ) = disable!(sampler, (:fabricate, 0))
#
# The fabrication rate depends on the number of parts present as a total,
# but each part can wear out individually with a Gamma distribution that starts
# its clock when that particular part is produced. The Gamma distribution is the
# same distribution we use to model time-to-wear-out.
#
function fabricate_part(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.parts)
    total = pre_event_total
    part_id = length(model.parts) + 1
    push!(model.parts, (when, true))
    time_offset = when - 0  # The 0 is when the machine turned on.
    rate1 = SetupRate(θ[:fabricate_max][1], θ[:fabricate_setup][1]; t0=time_offset)
    rate2 = SetupRate(θ[:fabricate_max][2], θ[:fabricate_setup][2]; t0=time_offset)
    enable!(sampler, (:fabricate, 0), [rate1, rate2])
    total = count(x -> x[2], model.parts)
    assemblerate1 = Exponential(inv(θ[:assemble][1] * total))
    assemblerate2 = Exponential(inv(θ[:assemble][2] * total))
    pre_event_total > 0 && disable!(sampler, (:assemble, 0))
    enable!(sampler, (:assemble, 0), [assemblerate1, assemblerate2])
    gamma1 = Gamma(θ[:scrap_k][1], inv(θ[:scrap_theta][1]))
    gamma2 = Gamma(θ[:scrap_k][2], inv(θ[:scrap_theta][2]))
    enable!(sampler, (:scrap, part_id), [gamma1, gamma2])
end
#
# Scrapping a part will turn off widget production for that part,
# but if it turns off the _last_ part, then it turns off all assembly.
#
function scrap_part(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.parts)
    total = pre_event_total
    model.parts[individual] = (zero(Time), false)
    pre_event_total > 0 && disable!(sampler, (:assemble, 0))
    total = count(x -> x[2], model.parts)
    if total > 0
        assemblerate1 = Exponential(inv(θ[:assemble][1] * total))
        assemblerate2 = Exponential(inv(θ[:assemble][2] * total))
        enable!(sampler, (:assemble, 0), [assemblerate1, assemblerate2])
    end
end
#
# When assembly happens, it isn't a one-time event. The same event
# can happen again, so we re-enable it.
#
function assemble_widget(model, sampler, individual, when, θ)
    pre_event_total = count(x -> x[2], model.parts)
    model.widgets += 1
    if pre_event_total > 0
        assemblerate1 = Exponential(inv(θ[:assemble][1] * pre_event_total))
        assemblerate2 = Exponential(inv(θ[:assemble][2] * pre_event_total))
        enable!(sampler, (:assemble, 0), [assemblerate1, assemblerate2])
    end
end
#
# Our central dynamics calls the appropriate function above when its event
# fires.
#
function step_workstation!(model, sampler, which, when)
    event, individual = which
    Dict(
        :on => start_machine,
        :off => stop_machine,
        :fabricate => fabricate_part,
        :scrap => scrap_part,
        :assemble => assemble_widget,
    )[event](model, sampler, individual, when, model.θ)
end
#
# The CompetingClocks.jl package does the sampling, but you create the mainloop.
#
function one_epoch(model, sampler)
    step_workstation!(model, sampler, (:on, 0), time(sampler))
    when, which = next(sampler)
    while !isnothing(which)
        fire!(sampler, which, when)
        step_workstation!(model, sampler, which, when)
        when, which = next(sampler)
    end
    basal, weighted = pathloglikelihood(sampler, time(sampler))
    logimportance = basal - weighted
    return (model.widgets, logimportance)
end
#
# This runs the simulation many times with chosen parameters.
#
function run_epochs(epoch_cnt, use_importance, rng)
    ## We define two sets of parameters. The first biases the simulation towards
    ## producing a rare event and the second is the basal rate we use to evaluate
    ## the importance of those events.
    params = Dict(
        :machine_off => (0.6, 0.2), # per minute
        :fabricate_max => (10.0, 10.0), # parts/min
        :fabricate_setup => (1.0, 1.0), # per minute, rate of machine warm-up.
        :scrap_k => (4.0, 4.0),  # k for Gamma
        :scrap_theta => (4 * 4 / 30, 4 * 4 / 30), # theta for Gamma
        :assemble => (1.00, 1.05), # widgets/min/part
    )
    widgets = zeros(Int, epoch_cnt)
    importance = zeros(Float64, epoch_cnt)
    model = Workstation(params)
    builder = SamplerBuilder(
        Tuple{Symbol,Int}, Float64;
        method=FirstToFireMethod(),
        path_likelihood=true,
        likelihood_cnt=2,
    )
    sampler = SamplingContext(builder, rng)
    sample_from_distribution!(sampler, use_importance ? 2 : 1)
    for epoch_idx in eachindex(widgets)
        (cnt, weight) = one_epoch(model, sampler)
        widgets[epoch_idx] = cnt
        importance[epoch_idx] = weight
        empty!(model)
        reset!(sampler)
    end
    return widgets, importance
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
# * This example follows a classic bursty-production model, where a machine
#   switches on and off and produces parts in intermittent bursts.
