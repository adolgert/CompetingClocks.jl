# # Transitions with Memory
#
# Some transitions can be paused and restarted. These transitions are said to have memory [Zimmerman:2007]. For instance, let's say there is an industrial process that takes a Gamma-distributed amount of time.
using Distributions
using Plots
bottler = Gamma(7.5, 1.0)
x = 0.0:0.05:12.0
plot(x, pdf.(bottler, x), xlims=(0.0, 12.0))
#
# Then that process is interrupted, maybe because the machinery is required for a different process with higher priority, and then we restart the industrial process. What is the distribution of the restarted process, assuming it picks up where it left off? It's shifted to the left. Shifting the distribution doesn't just translate the pdf. It rescales the pdf so that the area under the curve is one and also changes the shape.
#
survival = ccdf
hazard(dist, x) = pdf(dist, x) / survival(dist, x)
conditional_survival(dist, intermediate, x) = survival(dist, x) / survival(dist, intermediate)
shiftpdf(dist, intermediate, x) = hazard(dist, x + intermediate) * conditional_survival(dist, intermediate, x + intermediate)
remembered_age = 5.1
plot(x, shiftpdf.(bottler, remembered_age, x), xlims=(0.0, 12.0))
#
# Another way to think about transitions with memory is that the zero-time of the process is now in the past. For Fleck, that zero-time is called the *enabling time.* When we enable a transition, we give the enabling time as an absolute time in the simulation.
#
# ```
# enablingtime = currenttime - remembered_age
# enable!(sampler, clock, distribution, enablingtime, currenttime, RNG)`
# ```
#
# If a simulation wants to include transitions that have memory, that simulation needs to store the total time a transition has been enabled and then use that to set the enabling time when it restarts a transtition.
#
# ## References
#
# [Zimmerman:2007] Zimmermann, Armin. Stochastic discrete event systems. Springer, Berlin Heidelberg New York, 2007.
#
