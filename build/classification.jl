#   Unit classification.jl, part of PosDefManifold Package for julia language
#   v 0.1.3 - last update 7th of Mai 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DESCRIPTION
#   This Unit implements classification methods and related functions
#   useful in relation to the Riemannian Geometry on the manifold
#   of Symmetric Positive Definite (SPD) or Hermitian matrices.
#
#   CONTENT
#   1. Probability
# __________________________________________________________________

#  ------------------------
## 1. Probability
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
