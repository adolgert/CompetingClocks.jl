using Distributions
using Random
import Distributions: params, partype, mean, median, mode, var, skewness, kurtosis
import Distributions: pdf, logpdf, cdf, ccdf, quantile, mgf, cf
import Base: rand
import Random: rand!

struct Never{T<:Real} <: ContinuousUnivariateDistribution
    Never{T}() where {T <: Real} = new{T}()
end

Never() = Never{Float64}()

params(d::Never) = ()
partype(d::Never{T}) where {T<:Real} = T
mean(d::Never) = Inf
median(d::Never) = Inf
mode(d::Never) = Inf
var(d::Never) = Inf
skewness(d::Never{T}) where {T<:Real} = zero(T)
kurtosis(d::Never{T}) where {T<:Real} = zero(T)
pdf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
logpdf(d::Never{T}, x::Real) where {T<:Real} = -Inf
cdf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
ccdf(d::Never{T}, x::Real) where {T<:Real} = one(x)
quantile(d::Never{T}, q::Real) where {T<:Real} = Inf
mgf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
cf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
rand(rng::Random.AbstractRNG, d::Never) = Inf

function rand!(rng::Random.AbstractRNG, d::Never, arr::AbstractArray)
    arr .= Inf
end
