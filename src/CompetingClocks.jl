module CompetingClocks

const ContinuousTime = AbstractFloat

include("setofsets.jl")
include("prefixsearch/binarytreeprefixsearch.jl")
include("prefixsearch/cumsumprefixsearch.jl")
include("prefixsearch/keyedprefixsearch.jl")
include("distributions/lefttrunc.jl")
include("distributions/hazard.jl")
include("distributions/primal.jl")
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
# `export` statements to individual source files.
#
# The BLESSED public surface is the high-level SamplingContext layer plus the
# SamplerBuilder/spec machinery. Only those names are exported. The low-level
# sampler layer remains fully usable via qualified `CompetingClocks.X` access
# and is declared `public` (Julia >= 1.11) so tooling recognizes it as part of
# the supported API without polluting the exported namespace.
# ---------------------------------------------------------------------------

# --- Blessed, exported surface --------------------------------------------

# Context layer (the mainstream high-level entry point)
export SamplingContext, enable!, disable!, fire!, next, next_delayed, reset!
export isenabled, enabled, steploglikelihood, pathloglikelihood
export sample_from_distribution!, freeze_crn!, reset_crn!
export enabled_history, disabled_history, clone, copy_clocks!, keytype, timetype

# Builder & sampler specs (user-facing way to choose a sampler)
export SamplerBuilder, add_group!, available_samplers
export NextReactionMethod, DirectMethod, FirstReactionMethod, FirstToFireMethod
export RejectionMethod, PartialPropensityMethod, PetriMethod

# Delayed reactions
export Delayed

# Distributions
export Never

# --- Low-level developer surface (public, not exported) -------------------
# These names are the developer/framework-author layer. They are reachable as
# `CompetingClocks.<Name>` and declared `public` on Julia >= 1.11. The bare
# `public` keyword only parses on 1.11+, so it is VERSION-gated for the
# julia = "^1.10" compat bound.
const _PUBLIC_NAMES = (
    :SSA, :MultiSampler, :SamplerChoice, :choose_sampler,
    :CombinedNextReaction, :FirstReaction, :FirstToFire,
    :DirectCall, :DirectCallExplicit, :MultipleDirect,
    :PSSACR, :RSSA, :Petri,
    :sampling_space, :set_bound!, :set_global_bound_factor!,
    :TrackWatcher, :DebugWatcher, :TrajectoryWatcher, :MemorySampler,
    :PathLikelihoods, :absolute_enabling,
    :CommonRandom, :misses, :misscount,
    :DelayedState, :SetOfSets, :build_sampler,
    :split!, :jitter!,
)

@static if VERSION >= v"1.11"
    eval(Meta.parse("public " * join(_PUBLIC_NAMES, ", ")))
end

end
