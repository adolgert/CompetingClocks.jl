using Distributions: UnivariateDistribution
import Distributions

"""
    primal_distribution(d::UnivariateDistribution)

The boundary through which the sampler receives a distribution. The sampler only
ever draws numbers, so it must run on plain `Float64` internals; it should never
see an automatic-differentiation tracer type in a parameter.

For ordinary number-typed parameters this is the identity: the distribution is
returned unchanged. When a package extension is loaded (see
`ext/CompetingClocksForwardDiffExt.jl`), it strips AD tracer types (e.g.
`ForwardDiff.Dual`) from the parameters and rebuilds a primal distribution, so
the sampler keeps working while the differentiable copy stays on the likelihood
watcher.

Dispatch is two-level: `primal_distribution` computes the parameter element type
with `Distributions.partype` and forwards to `_primal_distribution`, which an
extension can specialize on that element type. Dispatching on `partype` (rather
than the distribution type) lets the extension add a method without overwriting
this identity fallback for ordinary number types.
"""
primal_distribution(d::UnivariateDistribution) = _primal_distribution(d, Distributions.partype(d))

# Identity for ordinary number-typed parameters. An extension specializes the
# second argument (the parameter element type) for AD tracer types.
_primal_distribution(d::UnivariateDistribution, ::Type) = d
