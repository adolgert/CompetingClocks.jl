module CompetingClocks

const ContinuousTime = AbstractFloat

include("setofsets.jl")
include("prefixsearch/binarytreeprefixsearch.jl")
include("prefixsearch/cumsumprefixsearch.jl")
include("prefixsearch/keyedprefixsearch.jl")
include("distributions/lefttrunc.jl")
include("distributions/hazard.jl")
include("sample/interface.jl")
include("sample/sampler.jl")
include("distributions/neverdist.jl")
include("trace/track.jl")
include("trace/trajectory.jl")
include("trace/path_likelihoods.jl")
include("trace/debug.jl")
include("sample/nrtransition.jl")
include("sample/firstreaction.jl")
include("sample/firsttofire.jl")
# include("sample/fixeddirect.jl")
include("sample/direct.jl")
include("sample/multiple_direct.jl")
include("sample/combinednr.jl")
include("sample/pssa_cr.jl")
include("sample/rssa.jl")
include("variance/with_common_random.jl")
include("sample/petri.jl")
include("samplerspec.jl")
include("delayed_state.jl")
include("sampler_builder.jl")
include("context.jl")

# ---------------------------------------------------------------------------
# Consolidated public API. ALL exports for the package live here; do not add
# `export` statements to individual source files. When you add a new public
# name, add it to the appropriate group below.
# ---------------------------------------------------------------------------

# Context layer (the mainstream high-level entry point)
export SamplingContext, isenabled, freeze_crn!, sample_from_distribution!
export enabled_history, disabled_history, next_delayed, timetype

# Builder & sampler specs (user-facing way to choose a sampler)
export SamplerBuilder, add_group!, build_sampler
export available_samplers
export NextReactionMethod, DirectMethod, FirstReactionMethod, FirstToFireMethod
export RejectionMethod, PartialPropensityMethod, PetriMethod

# Low-level samplers and their shared operations
export SSA, MultiSampler, SamplerChoice, choose_sampler
export CombinedNextReaction, FirstReaction, FirstToFire
export DirectCall, DirectCallExplicit, MultipleDirect
export PSSACR, RSSA, Petri
export enable!, disable!, fire!, next, enabled, clone, copy_clocks!, reset!
export sampling_space, set_bound!, set_global_bound_factor!

# Watchers & likelihood (tracing, debugging, path likelihoods)
export TrackWatcher, DebugWatcher, TrajectoryWatcher, MemorySampler
export PathLikelihoods, steploglikelihood, pathloglikelihood, absolute_enabling

# Common random numbers (variance reduction)
export CommonRandom, misses, misscount, reset_crn!

# Delayed reactions
export Delayed, DelayedState

# Utilities
export Never, SetOfSets
export getindex, keys, length, keytype

end
