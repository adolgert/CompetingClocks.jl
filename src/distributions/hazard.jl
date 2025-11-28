
function hazard(d::UnivariateDistribution, te::Float64, t::Float64)
    (t < te) && return zero(Float64)
    df = logpdf(d, t - te) - logccdf(d, t - te)
    return exp(df)
end


hazard(d::Exponential, te::Float64, t::Float64) = (t >= te) ? rate(d) : zero(Float64)
