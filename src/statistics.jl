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
    mean scalar function
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
