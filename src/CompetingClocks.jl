module CompetingClocks
using Documenter

const ContinuousTime = AbstractFloat

include("setofsets.jl")
include("prefixsearch/binarytreeprefixsearch.jl")
include("prefixsearch/cumsumprefixsearch.jl")
include("prefixsearch/keyedprefixsearch.jl")
include("lefttrunc.jl")
include("sample/interface.jl")
include("sample/sampler.jl")
include("sample/neverdist.jl")
include("trace/track.jl")
include("trace/trajectory.jl")
include("trace/debug.jl")
include("sample/nrtransition.jl")
include("sample/firstreaction.jl")
include("sample/firsttofire.jl")
# include("sample/fixeddirect.jl")
include("sample/direct.jl")
include("sample/multiple_direct.jl")
include("sample/combinednr.jl")
include("variance/with_common_random.jl")
include("sample/petri.jl")
include("context.jl")

export Never

end
