using Test
using SafeTestsets

using CompetingClocks
using Distributions
using Random

include("travel.jl")
include("generate_data.jl")
include("single_clock.jl")
include("experiments.jl")

using .TravelModel
