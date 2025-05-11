#   Unit statistics.jl, part of PosDefManifold Package for julia language
#
#   MIT License
#   Copyright (c) 2019-25, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DESCRIPTION
#   This Unit implements routines for statistics and probability.
#
#   CONTENT
#   1.  Utilities
#   2.  Probability
#   3.  Descriptive Statistics
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

 **Examples**
```julia
χ=[1.0, 2.3, 0.4, 5.0]
π=softmax(χ)
```

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

 **Examples**
```julia
using PosDefManifold
# Generate 10 random numbers distributed as a chi-square with 2 df.
ν=[randχ²(2) for i=1:10]
arithmetic_mean=mean(Euclidean, ν)
geometric_mean=mean(Fisher, ν)
harmonic_mean=mean(invEuclidean, ν)
harmonic_mean<=geometric_mean<=arithmetic_mean # AGH inequality
```

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

"""
    std(metric::Metric, ν::Vector{T};
        corrected::Bool=true,
        mean=nothing) where T<:RealOrComplex

 Standard deviation of ``k`` real or complex scalars,
 using the specified `metric`
 of type [Metric::Enumerated type](@ref) and the
 specified `mean` if provided.

 Only the Euclidean and Fisher
 metric are supported by this function. Using the Euclidean
 metric return the output of standard Julia
 [std](https://docs.julialang.org/en/v1/stdlib/Statistics/#Statistics.std)
 function. Using the Fisher metric return the scalar geometric standard deviation,
 which is defined such as,

 ``\\sigma=\\text{exp}\\Big(\\sqrt{k^{-1}\\sum_{i=1}^{k}\\text{ln}^2(v_i/\\mu})\\Big)``.

 If `corrected` is `true`, then the sum is scaled with ``k-1``,
 whereas if it is `false` the sum is scaled with ``k``.

 **Examples**
```julia
using PosDefManifold
# Generate 10 random numbers distributed as a chi-square with 2 df.
ν=[randχ²(2) for i=1:10]
arithmetic_sd=std(Euclidean, ν) # mean not provided
geometric_mean=mean(Fisher, ν)
geometric_sd=std(Fisher, ν, mean=geometric_mean) # mean provided
```

"""
function std(metric::Metric, ν::Vector{T};
             corrected::Bool=true,
             mean=nothing) where T<:RealOrComplex

    if      metric == Euclidean     return std(ν; corrected=corrected, mean=mean)

    elseif  metric == Fisher
            mean===nothing ? μ=mean(Fisher, ν) : μ=mean
            if corrected return exp(√(𝚺(log(w/μ)^2 for w in ν)/(length(ν)-1)))
            else         return exp(√(𝛍(log(w/μ)^2 for w in ν))) end

    elseif  metric in (invEuclidean, logEuclidean, Jeffrey,
                        logdet0, ChoEuclidean, logCholesky,
                        VonNeumann, Wasserstein)
        @warn "function statistics.mean (scalar mean) not implemented for metric $metric"
    else
        @error "in RiemannianGeometry.mean function: the chosen 'metric' does not exist"
    end # if metric
end
