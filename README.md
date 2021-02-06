# Fleck

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://adolgert.github.io/Fleck.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://adolgert.github.io/Fleck.jl/dev)
[![Build Status](https://github.com/adolgert/Fleck.jl/workflows/CI/badge.svg)](https://github.com/adolgert/Fleck.jl/actions)

This is a continuous-time simulation with the simplest-possible state and non-exponential transitions. The goal is to make a compact, fastest-possible non-exponential continuous-time simulation.

The simplest-possible state is a vector-addition system, which has a state that is a vector of integers just like chemical equations. A transition is defined by a vector of integers that represent how many values to take or give. However, for this simulation, the transition rates can be any function of the state and have exponential or non-exponential distributions.
