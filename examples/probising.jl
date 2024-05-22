# Probabilistic Ising Model
#
# There is a simple statistical model that informs our understanding of
# avalanches, earthquakes, and noisy non-equlillibrium systems. This
# model is a two-dimensional checkerboard, where each square on the board
# represents a magnet pointed either up or down. Each square exerts
# force on its neighbors to flip them, and the whole of the board
# is in an applied magnetic field, which pushes on all of the magnets.
# The squares tend to flip, one at a time, and influence their
# neighbors to flip. Each cascade of neighbors-flipping-neighbors
# is an avalanche.
#
# Here, we'll put the board into an increasing magnetic field
# and ask how the rate of change of the magnetic field affects
# the size of avalanches.
#
# The energy for a single magnet interacting with the applied
# field is -mH, where the H changes in time, so H=-H0+H1 t.
# The single magnet also sees its neighbors, which changes
# its energy by -m(sum over neighboring m_j). And each magnet
# has a slight internal preferential direction, m_0.
# 
# The hazard rate for a single magnet to flip is then proportional
# to
#
# \lambda = exp(kT - mH - m\sum_mj - mm_0)
#
# This hazard is time-dependent.

mutable struct Board
    moment::Matrix{Float64}
    orientation::Matrix{Int}
    neighbor_coupling::Float64
end


struct Experiment
    board::Board
    starting_field::Float64
    ramp_rate::Float64
end


function new_experiment(ramp_rate, coupling, sidelen, rng)
    
end
