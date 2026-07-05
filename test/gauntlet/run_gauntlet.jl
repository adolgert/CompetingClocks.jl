# run_gauntlet.jl - Interactive / batch driver for the milestone-1 gauntlet.
#
# This is NOT part of the default CI suite. It runs the gauntlet at a real
# sample size and prints the verdict table. Reproduce a run with, e.g.:
#
#   julia --project=test test/gauntlet/run_gauntlet.jl 1000 20260704
#
# where the first argument is n_replications and the second is the seed.

using CompetingClocks
using CompetingClocks: FirstReaction, SSA, TrackWatcher, NextReactionMethod, FirstReactionMethod
using Random
using Distributions
using Graphs
using Base.Threads
using HypothesisTests

include("travel.jl")
using .TravelModel
include("generate_data.jl")
include("mark_calibration.jl")
include("doob_meyer.jl")
include("ad_diagnostics.jl")
include("anderson_darling.jl")
include("running_score.jl")
include("experiments.jl")
include("runner.jl")

function main()
    n_replications = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1000
    seed = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20260704
    verdicts = run_gauntlet(
        NextReactionMethod();
        n_replications = n_replications,
        seed = seed,
        state_cnt = 5,
        history_steps = 5,
        verbose = true,
    )
    return verdicts
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
