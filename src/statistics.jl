#   Unit statistics.jl, part of PosDefManifold Package for julia language
#   v 0.3.1 - last update 30th of Mai 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DESCRIPTION
#   This Unit implements routines for statistics and probability.
#
#   CONTENT
#   1.  Utilities
#   2.  Probability
#   3.  Descriptive Statistocs
# __________________________________________________________________


#  ------------------------
## 1. Utilities
#  ------------------------


#  ------------------------
## 2. Probability
#  ------------------------

"""
    softmax(χ::Vector{T}) where T<:Real

 Given a real vector of ``k`` non-negative scores ``χ=c_1,...,c_k``,
 return the vector ``π=p_1,...,p_k`` of their
 [softmax](https://en.wikipedia.org/wiki/Softmax_function) probabilities,
 as per

 `` p_i=\\frac{\\textrm{e}^{c_i}}{\\sum_{i=1}^{k}\\textrm{e}^{c_i}} ``.

 ## Examples
    χ=[1.0, 2.3, 0.4, 5.0]
    π=softmax(χ)

"""
softmax(χ::Vector{T}) where T<:Real = exp.(χ) ./ 𝚺(exp.(χ))


#  -------------------------
## 3. Descriptive Statistics
#  -------------------------
"""
    mean(metric::Metric, ν::Vector{T}) where T<:RealOrComplex

 Mean of ``k`` real or complex scalars, using the specified `metric`
 of type [Metric::Enumerated type](@ref). Note that using the Fisher,
 logEuclidean and Jeffrey metric, the resulting mean
 is the scalar geometric mean. Note also that the code of this method
 is in unit *statistics.jl*, while the code for all the others is
 in unit *riemannianGeometry.jl*.

 ## Examples
    using PosDefManifold
    # Generate 10 random numbers distributed as a chi-square with 2 df.
    ν=[randχ²(2) for i=1:10]
    arithmeticMean=mean(Euclidean, ν)
    geometricMean=mean(Fisher, ν)
    HarmonicMean=mean(invEuclidean, ν)
    HarmonicMean<=geometricMean<=arithmeticMean # AGH inequality

"""
function mean(metric::Metric, ν::Vector{T}) where T<:RealOrComplex
    if      metric == Euclidean     return mean(ν)
    elseif  metric == invEuclidean  return inv(mean(inv, ν))
    elseif  metric in (Fisher,
                       logEuclidean,
                       Jeffrey)     return exp(mean(log, ν))
    elseif  metric == logdet0
        @warn "function statistics.mean (scalar mean) not implemented for metric $metric"
    elseif  metric == ChoEuclidean  return (mean(√, ν))^2
    elseif  metric == logCholesky   return exp((mean(log, sqrt.(ν))))^2
    elseif  metric == VonNeumann
        @warn "function statistics.mean and .geodesic not defined for metric $metric"
    elseif  metric == Wasserstein   return (mean(√, ν))^(-0.5)
    else
        @error "in RiemannianGeometry.mean function: the chosen 'metric' does not exist"
    end # if metric
end
