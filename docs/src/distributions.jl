# # Non-exponential Simulation

using DisplayAs #hide
using Plots #hide
using Distributions #hide
using LaTeXStrings #hide
using QuadGK #hide
using ForwardDiff #hide

# Fleck is a sampler for generalized semi-Markov processes (GSMP). Every event in a generalized semi-Markov process is chosen as the result of a competion among clocks to see which fires next.
# 
# In a *process-oriented* simulation (like [SimJulia](https://simjuliajl.readthedocs.io/en/stable/welcome.html)), control flow is based on tasks. Each task performs some action on the state, rolls the dice, and sets a wake-up time. It might wake up as expected and possible execute code, or it might be interrupted by another task's actions. In contrast, an *event-oriented* simulation using Fleck will create a set of possible next events, assign a probability distribution for *when* each can happen, and the timing of which happens first determines *which* next event happens. Let's look at how a probability distribution describes the time for an event to happen and then how they compete in Fleck.
# 
# ## Distributions in Time
# 
# Let's say you have a cold. You know you aren't going to recover immediately, but, as days go by, you're more and more sure you'll recover soon. You can read this graph below as the probability to recover given that you have not yet recovered.

const Gamma51 = Gamma(5.0, 1.0) #hide
survival = ccdf #hide
hazard(dist, x) = pdf(dist, x) / survival(dist, x) #hide
conditional_survival(dist, intermediate, x) = survival(dist, x) / survival(dist, intermediate) #hide
shiftpdf(dist, intermediate, x) = hazard(dist, x) * conditional_survival(dist, intermediate, x) #hide
x = 0:0.01:15 #hide
y = hazard.(Gamma51, x) #hide
p = plot(x, y, label="Gamma(5,1)") #hide
xlabel!(p, "Time") #hide
ylabel!(p, "Hazard Rate [per unit time]") #hide
title!(p, "Hazard of a Gamma Distribution") #hide
DisplayAs.PNG(p) #hide

# 
# It starts at zero, meaning there's no way you'll recover as soon as you're sick. It gets more likely over time that you're at the tail end of being sick. This is called a *hazard rate,* and the hazard rate shown is that of a Gamma distribution, commonly used to describe the rate of recovery for a population of individuals who are sick.
# 
# If, instead, you want to see the number of people who recover on any given day, that is called a probability distribution function.

x = 0:0.01:15 #hide
y = pdf.(Gamma51, x) #hide
p = plot(x, y, label="Gamma(5,1)") #hide
xlabel!(p, "Time") #hide
ylabel!(p, "Probability Density Function") #hide
title!(p, "PDF of a Gamma Distribution") #hide
DisplayAs.PNG(p) #hide

# 
# Where the hazard rate is an instantaneous quantity at a point in time, the probability distribution function (pdf) integrates over all possible future times. If we call the hazard rate $\lambda(t)$ and call the pdf $f$, we get this relationship.
# 
# ```math
# f(t) = \lambda(t) e^{-\int_0^t \lambda(s)ds}
# ```
# 
# The graph of the pdf tells us that the most likely time for this event is a little before time 5, in whatever units. You will see graphs of [pdfs on Wikipedia](https://en.wikipedia.org/wiki/Gamma_distribution#/media/File:Gamma_distribution_pdf.svg) because this is how people usually think about the probability an event happens at some time.
# 
# A simulation, however, has multiple events possible at any one time. One event may happen, and then other events need to restart. Let's ask, if you still have a cold on day 5, what is the probability distribution function for when you will recover?

t=5.0 #hide
x1 = 5:0.01:15 #hide
y1 = shiftpdf.(Gamma51, t, x1) #hide
sgp = plot(x1, y1, label="Shifted Gamma(5,1)") #hide
x = 0:0.01:15 #hide
y = pdf.(Gamma51, x) #hide
plot!(sgp, x, y, label="Gamma(5,1)") #hide
vline!(sgp, [t], label="New start time") #hide
xlabel!(sgp, "Time") #hide
ylabel!(sgp, "Probability Density Function") #hide
title!(sgp, "Change to PDF After Not Firing") #hide
DisplayAs.PNG(sgp) #hide

# 
# The probability distribution function changes now that you know you didn't recover earlier than day 5. On the other hand, the hazard rate for recovery from the cold will be unchanged. Using the same hazard rate, we can recalculate the pdf from the new time $t_0=5$.
# 
# ```math
# f(t;t>t_0) = \lambda(t) e^{-\int_{t_0}^t \lambda(s)ds}
# ```
# 
# In some sense, the underlying nature of when a particular event will happen is described more by the hazard rate, whereas the distribution function tells us about ensembles of events.
# 
# The hazard rate is related to the well known cumulative distribution function (CDF) by an integral. The CDF tells us what is the overall probability the event occured some time in the interval ``[t_0,t_1)``.
# 
# ```math
# F(t_0;t_1) = 1 - e^{-\int_{t_0}^{t_1} \lambda(s) ds}
# ```
# 
# Equally important for simulation is the survival function (sometimes called the complementary cumulative distribution function), which is the probability the event will not occur until after ``t_1``.
# 
# ```math
# S(t_1) = 1 - F(t_1;0) = e^{-\int_{0}^{t_1} \lambda(s) ds}
# ```

y1 = cdf.(Gamma51, x) #hide
y2 = survival.(Gamma51, x) #hide
p = plot(x, y1, label="CDF") #hide
plot!(p, x, y2, label="Survival") #hide
xlabel!(p, "Time") #hide
ylabel!(p, "Probability") #hide
title!(p, "Cumulative distribution and Survival functions \nof a Gamma Distribution") #hide
DisplayAs.PNG(p) #hide

# For our example, survival is the chance the cold lasts longer than the given time.

# ## Competition
# 
# ### Individual Distributions
# 
# Let's think of a moment when there are three possible next events. There is a Gamma distribution for when you recover from a cold, a Weibull distribution for when you decide to take medicine for the cold, and an Exponential distribution for when your Mom calls you. Each one is described by a distribution in time, and we can think of them as three hazard rates.
# 
com_dists = [ #hide
    Exponential(10), #hide
    Weibull(1.8, 8), #hide
    Gamma(5.0, 1.5) #hide
] #hide
com_diff_dists = [ #hide
    Exponential(10.0), #hide
    Weibull(1.8, 8), #hide
    Gamma{ForwardDiff.Dual}(5.0, 1.5) #hide
] #hide
com_labels = ["Exponential(10)" "Weibull(1.8, 8)" "Gamma(5,1.5)"] #hide
com_x = 0:0.01:15 #hide
com_hazard = stack([hazard.(d, com_x) for d in com_dists]) #hide
a = plot(com_x, com_hazard, #hide
    label=com_labels, #hide
    title="Hazards of Individual Distributions") #hide
com_pdf = stack([pdf.(d, com_x) for d in com_dists]) #hide
b = plot(com_x, com_pdf, #hide
    label=["Exponential(10)" "Weibull(1.8, 8)" "Gamma(5,1.5)"], #hide
    title="PDFs of Individual Distributions" #hide
    ) #hide
indp = plot(a, b, layout=(1, 2), size=(1000, 400)) #hide
DisplayAs.PNG(indp) #hide

# 
# The separate hazard rates are what we put into the simulation. Given their competition, the hazard rates will remain unchanged, but the pdfs will change.
# 
# ### Total Probability
# 
# Each of the three clock distributions above corresponds to a unique event ``E_i``, which has a probability that it will be the first to fire. We calculate this probability by marginalizing over the other events, which ends up being an integral over the distribution, multiplied by the survivals of the other events.
# 
# ```math
# P[E_i] = \int_0^\infty f_i(t) \prod_{j\ne i} S_j(t) dt
# ```
# 
# That gives the graph on the left.
# 
function marginal_probability_event(dists, which, t) #hide
    core_matrix = pdf(dists[which], t) #hide
    for i in eachindex(dists) #hide
        if i != which #hide
            core_matrix *= survival(dists[i], t) #hide
        end #hide
    end #hide
    core_matrix #hide
end #hide
#hide
total_probability(which, lim) = quadgk(t -> marginal_probability_event(com_dists, which, t), 0, lim, rtol=1e-3)[1] #hide
tp = total_probability.([1, 2, 3], Inf) #hide
@show tp #hide
bb = bar(reshape(com_labels, 3), tp, legend=:none, ylabel="Probability of Each Clock") #hide
com_pdfs = zeros(Float64, length(com_x), length(com_dists)) #hide
for pdfidx in eachindex(com_dists) #hide
    for (ix, x) in enumerate(com_x) #hide 
        com_pdfs[ix, pdfidx] = marginal_probability_event(com_dists, pdfidx, x) #hide
    end #hide
end #hide
bp = plot(com_x, com_pdfs, labels=com_labels, title="PDFs of Competing Clocks") #hide
condp = plot(bb, bp, size=(1000, 400)) #hide
DisplayAs.PNG(condp) #hide
# 
# The graph on the right shows the conditional distribution in time for each event, given that it was the one that fired, so it is $P[t_i | E_i]$. You can see that these distributions don't match the distributions for the individual events. They are modified by competition.
# 
# 
# ### Total Time
# 
# From the graph above, if we pick a time, $t_1=10$, we can read from the graph three hazard rates, $(\lambda_1(t_1),\lambda_2(t_1), \lambda_3(t_1))$. Each hazard rate is the rate, per unit time, of that event. We know that, if the simulation makes it to $t=10$ without any event happening, the conditional probability for any one of those events is the ratio of hazard rates.
# 
# ```math
# P[E_1|t=t_1] = \frac{\lambda_1(t_1)}{\lambda_1(t_1)+\lambda_2(t_1)+\lambda_3(t_1)}
# ```
# 
# If all distributions are Exponential, then they all have a constant hazard rate, and that conditional probability stays the same over time. If any distributions are non-Exponential, then that conditional probability changes.
# 
function total_pdf(t) #hide
    factor1 = zero(Float64) #hide
    for dist in com_dists #hide
        factor1 += hazard(dist, t) #hide
    end #hide
    factor2 = one(Float64) #hide
    for dist in com_dists #hide
        factor2 *= survival(dist, t) #hide
    end #hide
    return factor1 * factor2 #hide
end #hide
a = plot(com_x, total_pdf.(com_x), legend=:none, title="PDF for First to Fire") #hide
#hide
function com_cond_prob(i, t) #hide
    hazard(com_dists[i], t) / sum(hazard.(com_dists, t)) #hide
end #hide
com_rel_prob = stack([com_cond_prob.(i, com_x) for i in eachindex(com_dists)]) #hide
b = plot(com_x, com_rel_prob, labels=com_labels, ylabel="Conditional Probability", #hide
    title="Which Event Fires Given Firing Time") #hide
DisplayAs.PNG(plot(a, b, size=(1000, 400))) #hide
# 
# On the left of this graph is the pdf for the first event of the three to fire. We can see this as a marginal $P[t]$ and then the right-hand graph as the conditional $P[E_i|t]$.
# 
# ## Specification of a Simulation
# 
# We simulated using competing clocks when we measure using survival analysis. One is the inverse of the other.
# 
# A classic example of survival analysis is a hospital study where participants enter a trial for a drug, and doctors observe the date at which each patient recovers. How would you draw a probability distribution function for such a trial? You would draw a histogram where each week shows a count of how many recovered. What about the hazard rate for such a trial? Here, you would look only at the number of participants remaining in the trial and ask what fraction recover each week. That's a rate of recovery given that participants have not yet recovered.
# 
# There are many disciplines that use survival analysis to measure rates of processes. Epidemiology measures rates of contagion, death, and recovery. Demography measures time to have children or move to another country. Reliability studies the time for a part to break and be fixed. Drug discovery measures the time to approve a chemical for the next trial stage. For traffic flow, it's the time to the next intersection. For package delivery, it's the time to the next package scan. All of these do use, or can use, survival analysis to observe a system to simulate.
