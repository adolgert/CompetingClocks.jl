module CompetingClocksForwardDiffExt

# Extension that teaches the primal boundary how to strip ForwardDiff.Dual
# parameters out of a distribution. See `CompetingClocks.primal_distribution`.

using CompetingClocks
using Distributions
import ForwardDiff

# Strip a Dual to its value recursively so NESTED duals collapse to a plain
# number. ForwardDiff.hessian differentiates a gradient, so it hands us
# Dual{Tag2,Dual{Tag1,Float64}} and one `ForwardDiff.value` peels only one layer.
_value(x::ForwardDiff.Dual) = _value(ForwardDiff.value(x))
_value(x) = x

# Dispatch target for the primal boundary: reached only when the distribution's
# parameter element type is a Dual, so the sampler must not see this object.
function CompetingClocks._primal_distribution(d::UnivariateDistribution, ::Type{<:ForwardDiff.Dual})
    return _reconstruct_primal(d)
end

# Generic distribution families: rebuild the plain family from its params,
# stripping any dual parameters. The map leaves non-dual params (e.g. an integer
# count in Binomial) untouched so the family's own constructor still accepts them.
function _reconstruct_primal(d::UnivariateDistribution)
    wrapper = typeof(d).name.wrapper
    stripped = map(p -> p isa ForwardDiff.Dual ? _value(p) : p, Distributions.params(d))
    try
        return wrapper(stripped...)
    catch err
        throw(ArgumentError(
            "CompetingClocks could not build a primal (derivative-free) copy of a " *
            "$(nameof(typeof(d))) distribution to hand to the sampler. Pass a plain " *
            "distribution family whose keyword/positional constructor accepts " *
            "value-stripped parameters. Underlying error: $err"))
    end
end

# Truncated wraps an inner distribution; a params splat does not reconstruct the
# wrapper, so rebuild it explicitly, recursing on the inner distribution and
# stripping any dual bounds. `_value(nothing)` is `nothing`, so open bounds pass
# through unchanged.
function _reconstruct_primal(d::Distributions.Truncated)
    inner = CompetingClocks.primal_distribution(d.untruncated)
    return truncated(inner; lower=_value(d.lower), upper=_value(d.upper))
end

end # module
