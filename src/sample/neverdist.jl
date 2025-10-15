using Distributions
using Random
import Distributions: params, partype, mean, median, mode, var, skewness, kurtosis
import Distributions: pdf, logpdf, cdf, ccdf, quantile, mgf, cf
import Base: rand
import Random

struct Never{T<:Real} <: ContinuousUnivariateDistribution
    Never{T}() where {T<:Real} = new{T}()
end

Never() = Never{Float64}()

params(d::Never) = ()
partype(d::Never{T}) where {T<:Real} = T
mean(d::Never{T}) where {T} = typemax(T)
median(d::Never{T}) where {T} = typemax(T)
mode(d::Never{T}) where {T} = typemax(T)
var(d::Never{T}) where {T} = typemax(T)
skewness(d::Never{T}) where {T<:Real} = zero(T)
kurtosis(d::Never{T}) where {T<:Real} = zero(T)
pdf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
logpdf(d::Never{T}, x::Real) where {T<:Real} = typemin(T)
cdf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
ccdf(d::Never{T}, x::Real) where {T<:Real} = one(x)
quantile(d::Never{T}, q::Real) where {T<:Real} = typemax(T)
mgf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
cf(d::Never{T}, x::Real) where {T<:Real} = zero(x)
rand(rng::Random.AbstractRNG, d::Never{T}) where {T} = typemax(T)

function Random.rand!(rng::Random.AbstractRNG, d::Never{T}, arr::AbstractArray) where {T}
    arr .= typemax(T)
end
