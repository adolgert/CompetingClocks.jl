using InteractiveUtils

export available_samplers, PetriMethod
export NextReactionMethod, DirectMethod, FirstReactionMethod, FirstToFireMethod
export RejectionMethod, PartialPropensityMethod

# Define types a user can choose to specify which sampler they want.
# We make these auxiliary types so that the inner methods can be configurable
# without introducing lots of complexity.
abstract type SamplerSpec end

"""
    NextReactionMethod()

Uses Anderson's Modified Next Reaction method for distributions in an
Exponential class (Exponential, Weibull, Erlang) and the Next Reaction method
for other distributions. Yes, these two methods are mixed into one sampler because
Julia's holy traits pattern makes it efficient to use the fastest sampler on
a distribution-by-distribution basis. Because it reuses draws, this is the best choice if
you want to do variance reduction with Common Random Numbers.
"""
struct NextReactionMethod <: SamplerSpec end
(::NextReactionMethod)(K, T) = CombinedNextReaction{K,T}()


"""
    DirectMethod()
    DirectMethod(memory::Symbol, search::Symbol)

Use this to specify any Direct method for Exponential distributions. Defaults
to `memory=:remove` so it limits memory growth over time but `memory=:keep`
will be faster if the space of clock keys is limited. Defaults to
`search=:tree` for best performance for many enabled clocks but `search=:scan`
is faster for small numbers of clocks. The different kinds of methods, like
"Optimized Direct Methods" amount to using different computer science techniques
for scanning sums of hazard rates, and that's what the `search` algorithm lets
you choose.
"""
struct DirectMethod <: SamplerSpec
    memory_management::Symbol
    search_algorithm::Symbol
    DirectMethod() = new(:remove, :tree)
    function DirectMethod(symbols...)
        mem = nothing
        search = nothing
        for sym in symbols
            if sym ∈ (:keep, :remove)
                isnothing(mem) || error("Specify one of :keep or :remove to DirectMethod")
                mem = sym
            elseif sym ∈ (:tree, :array)
                isnothing(search) || error("Specify one of :tree or :array to DirectMethod")
                search = sym
            else
                error("Specify one of :keep, :remove, :tree, or :array to DirectMethod, not $sym")
            end
        end
        new(isnothing(mem) ? :remove : mem, isnothing(search) ? :tree : search)
    end
end
function (ss::DirectMethod)(K,T)
    DirectCallExplicit(
        K,
        T,
        ss.memory_management == :remove ? KeyedRemovalPrefixSearch : KeyedKeepPrefixSearch,
        ss.search_algorithm == :tree ? BinaryTreePrefixSearch : CumSumPrefixSearch,
        )
end


"""
    FirstReactionMethod()

The classic sampler that draws every clock at every time step. Very fast for
very small numbers of enabled clocks and returns a new `next()` every time it
is called which helps when resampling paths. Other samplers return the same
value every time you call `next()` unless you `jitter!` them, which is expensive.
"""
struct FirstReactionMethod <: SamplerSpec end
(::FirstReactionMethod)(K, T) = FirstReaction{K,T}()


"""
    FirstToFireMethod()

The simplest and fastest sampler. When you `enable!()` a clock, this draws the
firing time and saves it in a queue.
"""
struct FirstToFireMethod <: SamplerSpec end
(::FirstToFireMethod)(K, T) = FirstToFire{K,T}()


"""
    RejectionMethod(bound_factor=1.05)

A rejection-based algorithm. Only for exponential distributions, this may be
the fastest for large simulations. The `bound_factor>= 1.0` controls the default
upper bounds. Set to 1.0 for no rejections, which reduces this to the direct
method.
"""
struct RejectionMethod <: SamplerSpec
    bound_factor::Float64
    RejectionMethod(bound_factor=1.05) = new(bound_factor)
end
(rssa::RejectionMethod)(K, T) = RSSA{K,T}(; bound_factor=rssa.bound_factor)


"""
    PartialPropensityMethod()

Exact, continuous-time sampler using composition-rejection over groups.
Only for exponential distributions. This variant of partial-propensity
composition-rejection implements the sampler but not the reaction network. The
reaction-network can be implemented outside of the sampler.
"""
struct PartialPropensityMethod <: SamplerSpec
    ngroups::Int
    PartialPropensityMethod(ngroups=64) = new(ngroups)
end
(pssa::PartialPropensityMethod)(K, T) = PSSACR{K,T}(; ngroups=pssa.ngroups)


"""
    PetriMethod()
    PetriMethod(dt)

Samples by picking at random ignoring distributions. Good for testing rare
cases in simulations. Increments time `dt=1.0` by default. It's called "Petri"
because a Petri net model always chooses the next event at random.
"""
struct PetriMethod <: SamplerSpec
    dt::Float64
    PetriMethod() = new(1.0)
    PetriMethod(dt) = new(dt)
end
(ss::PetriMethod)(K, T) = Petri{K,T}(ss.dt)


"""
    available_samplers()

Returns a list of docstrings for all available samplers.
"""
function available_samplers()
    return [string(Base.Docs.doc(atype)) for atype in InteractiveUtils.subtypes(SamplerSpec)]
end
