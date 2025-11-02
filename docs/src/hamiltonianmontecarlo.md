# Hamiltonian Monte Carlo

There is an excellent paper that shows how to do Hamiltonian Monte Carlo using
the likelihood of a path of a simulation. The core idea is that, while you set
up the sequenece of events and define distributions for when next events happen,
you don't have to run the simulation. You can instead let an HMC sampler choose
times for events in the system.

With this package, that means the HMC updates a list of events, like:

```julia
events = [(evt=(:infect, 0, 1), time=0.1), (evt=(:infect, 0, 3), time=0.15)...] 
```

Then you pass that event list into a likelihood calculator.
```julia
function one_epoch(model, sampler, events_list)
    for idx in eachindex(events_list)
        which = events_list[idx].evt
        when = events_list[idx].time
        step_model!(model, sampler, which, when)
    end
    return pathloglikelihood(sampler, time(sampler))
end
```
Then use that likelihood to guide the next round of HMC. Using this simulation-plus-sampler
setup is much easier than working out the math of non-Exponential distributions.

 1. Billig, E. M., Roy, J. A., Ross, M. E., Dolgert, D. J., & Levy, M. Z. (2015, October). A BAYESIAN MODEL FOR IDENTIFYING AND PREDICTING THE DYNAMICS OF URBAN INSECT INFESTATIONS. In AMERICAN JOURNAL OF TROPICAL MEDICINE AND HYGIENE (Vol. 93, No. 4, pp. 537-537). 8000 WESTPARK DR, STE 130, MCLEAN, VA 22101 USA: AMER SOC TROP MED & HYGIENE.
