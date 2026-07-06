
# Times are typed `::Real` (not `::Float64`) so that ForwardDiff.Dual-
# parameterized distributions can flow through: the hazard's number type follows
# the distribution's parameter type, which is where a θ-gradient lives.
function hazard(d::UnivariateDistribution, te::Real, t::Real)
    # Return the zero of the value's number type so a `t < te` early return does
    # not narrow a Dual result back to Float64. It also avoids evaluating
    # `logpdf` at a negative argument, whose dual part would be NaN.
    R = promote_type(Distributions.partype(d), typeof(t - te))
    (t < te) && return zero(R)
    df = logpdf(d, t - te) - logccdf(d, t - te)
    return exp(df)
end


# Exponential has a constant hazard equal to its rate, so skip the
# logpdf/logccdf path. `zero(r)` (not `zero(Float64)`) keeps the result on the
# parameter's number type when the rate is a Dual.
function hazard(d::Exponential, te::Real, t::Real)
    r = rate(d)
    return (t >= te) ? r : zero(r)
end
