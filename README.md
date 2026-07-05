# CompetingClocks

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://adolgert.github.io/CompetingClocks.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://adolgert.github.io/CompetingClocks.jl/dev)
[![Build Status](https://github.com/adolgert/CompetingClocks.jl/workflows/CI/badge.svg)](https://github.com/adolgert/CompetingClocks.jl/actions)

CompetingClocks samples continuous-time probability distributions with time-varying hazard rates. It is a sampler for stochastic simulation in continuous time. You give it the probability distribution functions, and it tells you which fires next.

## Used by

 * [ChronoSim.jl](https://github.com/adolgert/ChronoSim.jl) — a simulation framework for discrete events in continuous time, built on this sampler.
 * [Quarton.jl](https://github.com/adolgert/Quarton.jl) — queueing theory models in continuous time.

This package was formerly called `Fleck`.
