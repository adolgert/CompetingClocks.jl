using DataStructures

struct MNRMTransition
    heap_handle::Int
    Lambda::Float64 # cumulative intensity
    distribution::UnivariateDistribution
    te::Float64  # Enabling time of distribution
    t0::Float64  # Enabling time of transition
end
