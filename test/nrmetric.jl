# This file supports the Combined Next Reaction sampler.
# It tests a distribution to determine whether it is better
# to sample it in linear space (with Next Reaction method) or
# log space (with Modified Next Reaction method). It judges
# based on two factors: accuracy and speed.

# I ran this with:
#
# cd test
# julia --project=. nrmetric.jl
# mv snippet*.tex ../docs/notes
# cd ../docs/notes
# latexmk -pdf distributions

using Distributions
using Printf
using StatsBase


function linear_accuracy!(buffer::AbstractArray, dist)
    for (i, v) in enumerate(LinRange(0.0001, 0.9999, length(buffer)))
        y = cquantile(dist, v)
        x = ccdf(dist, y)
        buffer[i] = abs(v - x)
    end
    sort!(buffer)
    average = mean(trim(buffer; prop=0.1))
    std = trimvar(buffer; prop=0.1)
    return (l1=buffer[end], mean=average, stdev=std)
end


function testing_minmax(dist)
    low = minimum(dist)
    if !isfinite(low)
        low = -10
    end
    high = maximum(dist)
    if !isfinite(high)
        high = 10
    end
    return (low, high)
end


function log_accuracy!(buffer::AbstractArray, dist)
    low, high = testing_minmax(dist)
    for (i, v) in enumerate(LinRange(low, high, length(buffer)))
        y = logccdf(dist, v)
        x = invlogccdf(dist, y)
        delta = abs(v - x)
        buffer[i] = delta
        if delta > 10e-6
            println("value $v logccdf $y inverse $x")
        end
    end
    sort!(buffer)
    average = mean(trim(buffer; prop=0.1))
    std = trimvar(buffer; prop=0.1)
    return (l1=buffer[end], mean=average, stdev=std)
end


function linear_timing!(buffer::AbstractArray, dist)
    # Use the inverse function first so that we generate times
    # that are in the support of the distribution.
    backward = @timed begin
        for (i, v) in enumerate(LinRange(0.0001, 0.9999, length(buffer)))
            buffer[i] = cquantile(dist, v)
        end
    end
    forward = @timed begin
        for i in eachindex(buffer)
            buffer[i] = ccdf(dist, buffer[i])
        end
    end
    return (forward, backward)
end


function log_timing!(buffer::AbstractArray, dist)
    # For the log-space, the forward-function is in support on [0, Inf)
    low, high = testing_minmax(dist)
    forward = @timed begin
        for (i, v) in enumerate(LinRange(low, high, length(buffer)))
            buffer[i] = logccdf(dist, v)
        end
    end
    backward= @timed begin
        for i in eachindex(buffer)
            buffer[i] = invlogccdf(dist, buffer[i])
        end
    end
    return (forward, backward)
end


function color(space, linear, logear, text)
    if !((space == "Linear") ⊻ (linear < logear))
        return raw"{\color{blue}" * text * raw"}"
    else
        return text
    end
end


"""
Test a distribution and print the results in a format for
a Latex table.
"""
function test_one!(buffer::AbstractArray, log::Vector{String}, dist)
    println("  linear")
    linear = (acc=linear_accuracy!(buffer, dist), perf=linear_timing!(buffer, dist))
    println("  log")
    logear = (acc=log_accuracy!(buffer, dist), perf=log_timing!(buffer, dist))

    name = repr(dist)
    for (unic, latexc) in [
        ("{", raw"\{"), ("}", raw"\}"),
        ("θ", raw"$\theta$"), ("α", raw"$\alpha$"), ("β", raw"$\beta$"),
        ("σ", raw"$\sigma$"), ("μ", raw"$\mu$", ("κ", raw"$\kappa$")),
        ("λ", raw"$\lambda$"), ("γ", raw"$\gamma$"), ("δ", raw"$\delta$"),
        ("χ", raw"$\chi$"), ("ξ", raw"$\xi$"), ("η", raw"$\eta$"),
        ("ν", raw"$\nu$")
        ]
        name = replace(name, unic => latexc)
    end
    push!(log, raw"\multicolumn{7}{|l|}{" * name * raw"}\\ \hline")
    for (space, metric) in [("Linear", linear), ("Log", logear)]
        l1 = @sprintf "%3.1f" log10(metric.acc.l1)
        eps = @sprintf "%3.1f" log10(metric.acc.mean)
        forward = @sprintf "%7.2e" metric.perf[1].time
        backward = @sprintf "%7.2e" metric.perf[2].time
        both = @sprintf "%7.2e" (metric.perf[1].time + metric.perf[2].time)
        l1 = color(space, linear.acc.l1, logear.acc.l1, l1)
        eps = color(space, linear.acc.mean, logear.acc.mean, eps)
        forward = color(space, linear.perf[1].time, logear.perf[1].time, forward)
        backward = color(space, linear.perf[2].time, logear.perf[2].time, backward)
        both = color(
            space,
            linear.perf[1].time + linear.perf[2].time,
            logear.perf[1].time + logear.perf[2].time,
            both
            )
        push!(log, "& $space & $l1 & $eps & $forward & $backward & $both" * raw"\\\\")
    end
    push!(log, raw"\hline")
end

sample_cnt = 100
buffer = zeros(Float64, sample_cnt)
test_distributions = [
    Distributions.Arcsine(),
    Distributions.BetaPrime(),
    Distributions.Biweight(),
    Distributions.Beta(),
    Distributions.Cauchy(),
    # Distributions.Chernoff(), # error newton not defined.
    # Distributions.Chi(), # no symbol
    # Distributions.Chisq(), # no symbol
    Distributions.Cosine(),
    Distributions.Epanechnikov(),
    Distributions.Erlang(),
    Distributions.Exponential(),
    # Distributions.FDist(),
    Distributions.Frechet(),
    Distributions.Gamma(),
    # Distributions.GeneralizedExtremeValue(), # no symbol
    Distributions.GeneralizedPareto(),
    Distributions.Gumbel(),
    Distributions.InverseGamma(),
    Distributions.InverseGaussian(),
    Distributions.JohnsonSU(),
    Distributions.Kolmogorov(),
    Distributions.Kumaraswamy(),
    Distributions.Laplace(),
    Distributions.Levy(),
    Distributions.Lindley(),
    Distributions.Logistic(),
    Distributions.LogitNormal(),
    Distributions.LogNormal(),
    # Distributions.LogUniform(), # no symbol
    # Distributions.NoncentralBeta(), # no symbol
    # Distributions.NoncentralChisq(), # no symbol
    # Distributions.NoncentralF(), # no symbol
    # Distributions.NoncentralT(), # no symbol
    Distributions.Normal(),
    Distributions.NormalCanon(),
    # Distributions.NormalInverseGaussian(),
    Distributions.Pareto(),
    Distributions.PGeneralizedGaussian(),
    Distributions.Rayleigh(),
    Distributions.Rician(),
    # Distributions.Semicircle(), # no symbol
    Distributions.SkewedExponentialPower(),
    # Distributions.SkewNormal(), # no method matching iterate(::SkewNormal{Float64})
    # Distributions.StudentizedRange(), # no symbol
    Distributions.SymTriangularDist(),
    # Distributions.TDist(), # no symbol
    # Distributions.TriangularDist(), # no symbol
    Distributions.Triweight(),
    Distributions.Uniform(),
    # Distributions.VonMises(), # no method matching iterate(::VonMises{Float64})
    Distributions.Weibull(),
]

println("There are $(length(test_distributions)) distributions to test.")

batch_size = 12
full, partial = divrem(length(test_distributions), batch_size)
iter_cnt = full + (partial > 0 ? 1 : 0)
for snip_idx in 1:iter_cnt
    log = Vector{String}()
    push!(log, raw"\begin{tabular}{|llrrrrr|} \hline")
    push!(log, raw"& Space & $\mbox{log}_{10}|\epsilon|_1$ & $\mbox{log}_{10}\langle\epsilon\rangle$ & forward [s] & backward [s] & both [s]\\ \hline")

    low = (snip_idx - 1) * batch_size + 1
    high = min(snip_idx * batch_size, length(test_distributions))
    for dist in test_distributions[low:high]
        println("Testing $(repr(dist))")
        test_one!(buffer, log, dist)
    end

    push!(log, raw"\end{tabular}")
    write("snippet$(snip_idx).tex", join(log, "\n"))
end
