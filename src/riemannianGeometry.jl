#    Unit riemannianGeometry.jl, part of PosDefManifold Package for julia language
#
#    MIT License
#    Copyright (c) 2019-25, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home
#
#    DESCRIPTION
#    This Unit implements Riemannian Geometry operations on the manifold
#    of Symmetric Positive Definite (SPD) or Hermitian matrices
#
#    CONTENT
#    0. Internal Functions
#    1. Geodesic Equations
#    2. Distances
#    3. Inter-distance matrices, Laplacian and Spectral Embedding
#    4. Means (center of mass, barycenters, ...)
#    5. Tangent Space
#    6. Procrustes Problems



# -----------------------------------------------------------
# 0. Internal Functions (not exported)
#    By convention their name begin with underscore char
# -----------------------------------------------------------

# Given a non-negative weight vector normalize the weights so as to sum up to 1
# if ✓w == true and if they are not already normalized
function _getWeights(w::Vector, ✓w::Bool)
    if ✓w==true
        s=𝚺(w)
        if s ≉  1.0 return w./s else return w end
    else return w
    end # if
end

# Set variables used by all iterative algorithms :
# k     = number of input matrices in `𝐏`
# n     = dimension of these matrices
# type  = type of the elements of these matrices
# thr   = number of threads set in Julia for multi-threading
# n²    = `n`²
# iter  = count of iterations. It is just initialized to 1
# conv  = convergence attained at each iteration. It is just initialized to 0
# oldconv = convergence attained at last iteration. It is just initiaized to maxpos=1e15
# converged = true if the algorithm has converged. It is just initialized to false
# ς     = step size for gradient descent. It is not used by all algorithms. It is just initialized to 1
# threaded = it is set to true if `k`>=`thr`*4 && `thr` > 1, false otherwise
# tolerance = tolerance for convergence. It depends on `type`. See help file of the concerned iterative algorithm
# v = vector of weights, either empty or initialized by function `_getWeights(w::Vector, ✓w::Bool)`
function _setVar_IterAlg(𝐏::Union{ℍVector, 𝔻Vector}, w::Vector, ✓w::Bool, tol::Real, ⏩::Bool)
    k, n, type, thr = dim(𝐏, 1), dim(𝐏, 2), eltype(𝐏[1]), nthreads()
    n², iter, conv, oldconv, converged, ς = n^2, 1, 0., maxpos, false, 1.
    ⏩ && thr > 1 ? threaded=true : threaded=false
    tol≈0. ? tolerance = √eps(real(type)) : tolerance = tol
    isempty(w) ? v=[] : v = _getWeights(w, ✓w)
    return (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v)
end

# used by all iterative algorithms.
# Print info to the user at the beginning of the algorithm if `verbose` is true
function _giveStartInfo_IterAlg(threaded::Bool, verbose::Bool, msg::String)
    if verbose
        println("")
        threaded && @info("Iterating multi-threaded "*msg*"...")
        !threaded && @info("Iterating "*msg*"...")
    end
end

# used by all iterative algorithms.
# print info to the user at the end of the algorithm if `verbose` is true
function _giveEndInfo_IterAlg(converged::Bool, verbose::Bool)
    verbose ? (converged ? @info("Convergence has been attained.\n") : @warn("Convergence has not been attained.")) : nothing
    verbose && println("")
end

# -----------------------------------------------------------
# 1. Geodesic Equations
# -----------------------------------------------------------

"""
```julia
    (1) geodesic(metric::Metric, P::ℍ{T}, Q::ℍ{T}, a::Real) where T<:RealOrComplex
    (2) geodesic(metric::Metric, D::𝔻{S}, E::𝔻{S}, a::Real) where S<:Real
```

(1) Move along the [geodesic](@ref) from point ``P`` to point ``Q``
(two positive definite matrices) with *arclegth* ``0<=a<=1``,
using the specified metric, of type [Metric::Enumerated type](@ref).

For all metrics,
- with ``a=0`` we stay at ``P``,
- with ``a=1`` we move up to ``Q``,
- with ``a=1/2`` we move to the mid-point of ``P`` and ``Q`` (mean).

Using the Fisher metric, argument ``a`` can be *any* real number, for instance:
- with ``0<a<1`` we move toward ``Q`` (*attraction*),
- with ``a>1`` we move over and beyond ``Q`` (*extrapolation*),
- with ``a<0`` we move back away from Q (*repulsion*).

``P`` and ``Q`` must be flagged by julia as `Hermitian`.
See [typecasting matrices](@ref).

The Fisher geodesic move is computed by the Cholesky-Schur algorithm
given in Eq. 4.2 by Iannazzo(2016)[🎓](@ref). If ``Q=I``,
the Fisher geodesic move is simply ``P^a`` (no need to call this funtion).

!!! note "Nota Bene"
    For the [logdet zero](@ref) and [Jeffrey](@ref) metric no closed form expression
    for the geodesic is available to the best of authors' knowledge,
    so in this case the geodesic is found as the weighted mean using the
    [`mean`](@ref) function.
    For the [Von Neumann](@ref) not even an expression for the mean is available,
    so in this case the geodesic is not provided and a *warning* is printed.

(2) Like in (1), but for two real positive definite diagonal matrices
``D`` and ``E``.

**Maths**

For points ``P``, ``Q`` and arclength ``a``, letting ``b=1-a``,
the geodesic equations for the supported metrics are:

| Metric   | geodesic equation |
|:----------:|:----------- |
|Euclidean| ``bP + aQ`` |
|invEuclidean| ``\\big(bP^{-1} + aQ^{-1}\\big)^{-1}``|
|ChoEuclidean| ``TT^*``, where ``T=bL_P + aL_Q``|
|logEuclidean| ``\\text{exp}\\big(b\\hspace{2pt}\\text{log}(P) + a\\hspace{2pt}\\text{log}(Q)\\big)``|
|logCholesky| ``TT^*``, where ``T=S_P+a(S_Q-S_P)+D_P\\hspace{2pt}\\text{exp}\\big(a(\\text{log}D_Q-\\text{log}D_P)\\big)``|
|Fisher | ``P^{1/2} \\big(P^{-1/2} Q P^{-1/2}\\big)^a P^{1/2}``|
|logdet0| uses weighted mean algorithm [`logdet0Mean`](@ref) |
|Jeffrey | uses weighted mean [`mean`](@ref) |
|VonNeumann | N.A.|
|Wasserstein| ``b^2P+a^2Q +ab\\big[(PQ)^{1/2} +(QP)^{1/2}\\big]``|

  **legend:** ``L_X``, ``S_X`` and ``D_X``
   are the Cholesky lower triangle of ``X``, its strictly lower triangular part
   and diagonal part, respectively (hence, ``S_X+D_X=L_X``,  ``L_XL_X^*=X``).

 **See also**: [`mean`](@ref).

 **Examples**
```julia
using PosDefManifold
P=randP(10)
Q=randP(10)
# Wasserstein mean
M=geodesic(Wasserstein, P, Q, 0.5)
# extrapolate suing the Fisher metric
E=geodesic(Fisher, P, Q, 2)
```
"""
function geodesic(metric::Metric, P::ℍ{T}, Q::ℍ{T}, a::Real) where T<:RealOrComplex
    if a ≈ 0 return P end
    if a ≈ 1 return Q end
    b = 1-a

    if      metric==Euclidean    return ℍ(P*b + Q*a)

    elseif  metric==invEuclidean return inv( ℍ(inv(P)b + inv(Q)a) )

    elseif  metric==logEuclidean return ℍ( exp( ℍ(log(P)b + log(Q)a) ) )

    elseif  metric==Fisher
            # Cholesky-Schur form (faster):
            L = cholesky(P, check=true)
            U⁻¹ = inv(L.U)
            F = schur(U⁻¹' * Q * U⁻¹)
            return ℍ(L.U' * (F.Z * F.T^a * F.Z') * L.U)

            # classical form (slower):
            #P½, P⁻½ = pow(P, 0.5, -0.5)
            #return ℍ( P½ * (P⁻½ * Q * P⁻½)^a * P½ )

    elseif metric ∈ (logdet0, Jeffrey)
            return mean(metric, ℍVector([P, Q]), w=[b, a], ✓w=false)

    elseif  metric==VonNeumann
            @warn("An expression for the geodesic is not available for the Von Neumann metric")

    elseif  metric==ChoEuclidean
            Z=choL(P)b + choL(Q)a
            return ℍ(Z*Z')

    elseif  metric==logCholesky
            LP=choL(P)
            LQ=choL(Q)
            slLP=tril(LP,-1)
            Z=slLP + a*(tril(LQ,-1)-slLP) +𝑓𝔻(x->x, LP)*exp((𝑓𝔻(log, LQ)a-𝑓𝔻(log, LP)a))
            return ℍ(Z*Z')

    elseif  metric==Wasserstein
            if T<:Real # isreal(P) && isreal(Q)
                    return ℍ( (b^2)*P + (a^2)*Q + (a*b)*real(√(P*Q)+√(Q*P)) )
            else    return ℍ( (b^2)*P + (a^2)*Q + (a*b)*(√(P*Q)+√(Q*P)) )
            end

    else    @error("in RiemannianGeometryP.geodesic function
                 (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end # function

function geodesic(metric::Metric, D::𝔻{T}, E::𝔻{T}, a::Real) where T<:Real
    if      a ≈ 0 return D end
    if      a ≈ 1 return E end
    b = 1-a

    if      metric==Euclidean    return D*b + E*a

    elseif  metric==invEuclidean return inv( inv(D)b + inv(E)a )

    elseif  metric ∈ (Fisher,
                 logEuclidean)   return exp( log(D)b + log(E)a )

    elseif  metric ∈ (logdet0,
                       Jeffrey)  return mean(metric, 𝔻Vector([D, E]), w=[b, a], ✓w=false)

    elseif  metric==VonNeumann
            @warn("An expression for the geodesic is not available for the Von Neumann metric")

    elseif  metric==ChoEuclidean
            Z=ℍ(√D)b + ℍ(√E)a;  return Z*Z

    elseif  metric==logCholesky # ???
            LD=ℍ(√D)
            Z=𝑓𝔻(x->x, LD)*exp((𝑓𝔻(log, ℍ(√E))a-𝑓𝔻(log, LD)a)); return Z*Z

    elseif  metric==Wasserstein  return (b^2)D + (a^2)E + (a*b)(D*E)

    else    @error("in RiemannianGeometryP.geodesic function
                 (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end # function



# -----------------------------------------------------------
# 2. Distances
# -----------------------------------------------------------

"""
```julia
    (1) distanceSqr(metric::Metric, P::ℍ{T}) where T<:RealOrComplex
    (2) distanceSqr(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    (3) distanceSqr(metric::Metric, D::𝔻{S}) where S<:Real
    (4) distanceSqr(metric::Metric, D::𝔻{S}, E::𝔻{S}) where S<:Real
```

**alias**: `distance²`

(1) Return ``δ^2(P, I)``, the *square of the distance* (or *divergence*) of positive definite
matrix ``P`` from the the identity matrix. See [distance from the origin](@ref).

(2) Return ``δ^2(P, Q)``, the *square of the distance* (or *divergence*) between two
positive definite matrices ``P`` and ``Q``. See [distance](@ref).

In both cases the distance function ``δ`` is induced by the argument `metric` of type
[Metric::Enumerated type](@ref).

``P`` in (1) and ``P``, ``Q`` in (2) must be flagged by julia as `Hermitian`.
See [typecasting matrices](@ref).

(3) and (4) are specialized methods of (1) and (2), respectively,
for real positive definite `Diagonal` matrices.
See [ℍVector type](@ref) and [𝔻Vector type](@ref).

**Maths**

For point ``P`` the *squared distances from the identity*
for the supported metrics are:

| Metric   | Squared Distance from the identity |
|:----------:|:----------- |
|Euclidean | ``∥P-I∥^2`` |
|invEuclidean| ``∥P^{-1}-I∥^2``|
|ChoEuclidean| ``∥L_P-I∥^2``|
|logEuclidean| ``∥\\textrm{log}P∥^2`` |
|logCholesky| ``∥S_P∥^2+∥\\textrm{log}D_P∥^2`` |
|Fisher | ``∥\\textrm{log}P∥^2`` |
|logdet0| ``\\textrm{logdet}\\frac{1}{2}(P+I) - \\frac{1}{2}\\textrm{logdet}(P)`` |
|Jeffrey | ``\\frac{1}{2}\\textrm{tr}(P+P^{-1})-n `` |
|VonNeumann | ``\\frac{1}{2}\\textrm{tr}(P\\textrm{log}P-\\textrm{log}P)``|
|Wasserstein| ``\\textrm{tr}(P+I) -2\\textrm{tr}(P^{1/2})`` |

For points ``P`` and ``Q`` their *squared distances* for the supported metrics are:

| Metric   | Squared Distance |
|:----------:|:----------- |
|Euclidean | ``∥P-Q∥^2`` |
|invEuclidean| ``∥P^{-1}-Q^{-1}∥^2``|
|ChoEuclidean| ``∥ L_P - L_Q ∥^2``|
|logEuclidean| ``∥\\textrm{log}P-\\textrm{log}Q∥^2`` |
|logCholesky| ``∥S_P-S_Q∥^2+∥\\textrm{log}D_P-\\textrm{log}D_Q∥^2`` |
|Fisher | ``∥\\textrm{log}(P^{-1/2}QP^{-1/2})∥^2``|
|logdet0| ``\\textrm{logdet}\\frac{1}{2}(P+Q) - \\frac{1}{2}\\textrm{logdet}(PQ)`` |
|Jeffrey | ``\\frac{1}{2}\\textrm{tr}(Q^{-1}P+P^{-1}Q)-n `` |
|VonNeumann | ``\\frac{1}{2}\\textrm{tr}(P\\textrm{log}P-P\\textrm{log}Q+Q\\textrm{log}Q-Q\\textrm{log}P)``|
|Wasserstein| ``\\textrm{tr}(P+Q) -2\\textrm{tr}(P^{1/2}QP^{1/2})^{1/2}`` | ``\\textrm{tr}(P+Q) -2\\textrm{tr}(P^{1/2}QP^{1/2})^{1/2}`` |

**legend:** ``L_X``, ``S_X`` and ``D_X``
are the Cholesky lower triangle of ``X``, its strictly lower triangular part
and diagonal part, respectively (hence, ``S_X+D_X=L_X``,  ``L_XL_X^*=X``).

**See also**: [`distanceSqrMat`](@ref).

**Examples (1)**
```julia
using PosDefManifold
P=randP(10)
d=distanceSqr(Wasserstein, P)
e=distanceSqr(Fisher, P)
metric=Metric(Int(logdet0)) # or metric=logdet0
s=string(metric) # check what is the current metric
f=distance²(metric, P) #using the alias distance²
```

**Examples (2)**
```julia
using PosDefManifold
P=randP(10)
Q=randP(10)
d=distanceSqr(logEuclidean, P, Q)
e=distance²(Jeffrey, P, Q)
```
"""
function distanceSqr(metric::Metric, P::ℍ{T}) where T<:RealOrComplex
    z=real(T)(0)

    if      metric==Euclidean    return max(z, ss(P-I))

    elseif  metric==invEuclidean return max(z, ss(inv(P)-I))

    elseif  metric ∈ (Fisher,
                logEuclidean)    return max(z, 𝚺(log.(eigvals(P)).^2))

    elseif  metric==logdet0      return max(z, real(logdet(0.5*(P+I)) - 0.5*logdet(P)))

    elseif  metric==ChoEuclidean return max(z, ss(choL(P)-I))

    elseif  metric==logCholesky
            LP=choL(P);          return max(z, real(sst(LP, -1)) + ssd(𝑓𝔻(log, LP)))

    elseif  metric==Jeffrey      return max(z, 0.5*(tr(P) + tr(inv(P))) - size(P, 1))

    elseif  metric==VonNeumann
            𝓵P=ℍ(log(P));        return max(z, 0.5*(tr(P, 𝓵P) - tr(𝓵P)))

    elseif  metric==Wasserstein  return max(z, tr(P) + size(P, 1) - 2*tr(sqrt(P)))

    else    @error("in RiemannianGeometryP.distanceSqr function
             (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end #function


function distanceSqr(metric::Metric, D::𝔻{T}) where T<:Real
    z=T(0)

    if     metric==Euclidean    return  max(z, ssd(D-I))

    elseif metric==invEuclidean return  max(z, ssd(inv(D)-I))

    elseif metric ∈ (Fisher,
                logEuclidean)   return  max(z, ssd(log(D)))

    elseif metric==logdet0      return  max(z, logdet(0.5*(D+I)) - 0.5*logdet(D))

    elseif metric==ChoEuclidean return  max(z, ssd(√D-I))

    elseif metric==logCholesky  return  max(z, ssd(𝑓𝔻(log, √D)))

    elseif metric==Jeffrey      return  max(z, 0.5*(tr(D) + tr(inv(D))) - size(D, 1))

    elseif metric==VonNeumann
           𝓵D=log(D);           return  max(z, 0.5*(tr(D*𝓵D) - tr(𝓵D)))

    elseif metric==Wasserstein  return  max(z, tr(D) + size(D, 1) - 2*tr(sqrt(D)))

    else   @error("in RiemannianGeometryP.distanceSqr function
             (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end #function


function distanceSqr(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    z=real(T)(0)
    if     metric==Euclidean    return  max(z, ss(ℍ(P - Q)))

    elseif metric==invEuclidean return  max(z, ss(ℍ(inv(P) - inv(Q))))

    elseif metric==logEuclidean return  max(z, ss(ℍ(log(P) - log(Q))))

    elseif metric==Fisher       return  max(z, 𝚺(log.(eigvals(P, Q)).^2))

    elseif metric==logdet0      return  max(z, real(logdet(0.5*(P + Q)) - 0.5*logdet(P * Q)))

    elseif metric==ChoEuclidean return  max(z, ss(choL(P)-choL(Q)))

    elseif metric==logCholesky
           LP=choL(P); LQ=choL(Q);
                                return  max(z, real(sst(tril(LP, -1) - tril(LQ, -1), -1)) + ssd(𝑓𝔻(log, LP) - 𝑓𝔻(log, LQ)))

    elseif metric==Jeffrey      return  max(z, 0.5*(tr(inv(Q), P) + tr(inv(P), Q)) - size(P, 1)) #using formula tr(Q⁻¹P)/2 + tr(P⁻¹Q)/2 -n

    elseif metric==VonNeumann # using formula: tr(PlogP - PlogQ + QlogQ - QlogP)/2=(tr(P(logP - LoqQ)) + tr(Q(logQ - logP)))/2=
           R=ℍ(log(P)-log(Q));  return  max(z, 0.5*real(tr(P, R) - tr(Q, R)))  # (tr(P(logP - LoqQ)) - tr(Q(logP - LoqQ)))/2

    elseif metric==Wasserstein
           P½=sqrt(P);          return  max(z, tr(P) + tr(Q) - 2*real(tr(sqrt(ℍ(P½ * Q * P½)))))

    else   @error("in RiemannianGeometryP.distanceSqr function
                    (PosDefManifold Package): the chosen 'metric' does not exist")
    end #if
end # function


function distanceSqr(metric::Metric, D::𝔻{T}, E::𝔻{T}) where T<:Real
    z=T(0)
    if     metric==Euclidean    return  max(z, ssd(D - E))

    elseif metric==invEuclidean return  max(z, ssd(inv(D) - inv(E)))

    elseif metric ∈ (Fisher,
                 logEuclidean)  return  max(z, ssd(log(D) - log(E)))

    elseif metric==logdet0      return  max(z, logdet(0.5*(D + E)) - 0.5*logdet(D * E))

    elseif metric==ChoEuclidean return  max(z, ssd(√(D) - √(E)))

    elseif metric==logCholesky  return  max(z, ssd(𝑓𝔻(log, √(D)) - 𝑓𝔻(log, √(E))))

    elseif metric==Jeffrey      return  max(z, 0.5*(tr(inv(E) * D) + tr(inv(D) * E)) - size(D, 1))

    elseif metric==VonNeumann
           R=log(D)-log(E);     return  max(z, 0.5*(tr(D * R) - tr(E * R)))

    elseif metric==Wasserstein  return  max(z, tr(D) + tr(E) - 2*tr(sqrt(D*E)))

    else   @error("in RiemannianGeometryP.distanceSqr function
                    (PosDefManifold Package): the chosen 'metric' does not exist")
    end #if
end # function
distance²=distanceSqr # alias


"""
```julia
    (1) distance(metric::Metric, P::ℍ{T}) where T<:RealOrComplex
    (2) distance(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    (3) distance(metric::Metric, D::𝔻{S}) where S<:Real
    (4) distance(metric::Metric, D::𝔻{S}, E::𝔻{S}) where S<:Real
```

(1) Return ``δ(P, I)``, the *distance* between positive definite matrix ``P`` and
the identity matrix.

(2) Return ``δ(P, Q)``, the *distance* between positive definite
matrices ``P`` and ``Q``.

(3) and (4) are specialized methods of (1) and (2), respectively,
for real positive definite `Diagonal` matrices.

This is the square root of [`distanceSqr`](@ref)
and is invoked with the same syntax therein.

**See also**: [`distanceMat`](@ref).
"""
distance(metric::Metric, P::ℍ{T}) where T<:RealOrComplex = √(distance²(metric, P))
distance(metric::Metric, D::𝔻{T}) where T<:Real = √(distance²(metric, D))

distance(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex = √(distance²(metric, P, Q))
distance(metric::Metric, D::𝔻{T}, E::𝔻{T}) where T<:Real = √(distance²(metric, D, E))



# -----------------------------------------------------------
# 3. Inter-distance matrix, Laplacian and Spectral Embedding
# -----------------------------------------------------------

"""
```julia
    (1) distanceSqrMat(metric::Metric, 𝐏::ℍVector;
    <⏩=true>)

    (2) distanceSqrMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;
    <⏩=true>) where T<:AbstractFloat
```

**alias**: `distance²Mat`

Given a 1d array ``𝐏`` of ``k`` positive definite matrices
``{P_1,...,P_k}`` of [ℍVector type](@ref), create the ``k⋅k`` real
`LowerTriangular` matrix comprising elements ``δ^2(P_i, P_j)\\textrm{, for all }i>=j``.

This is the lower triangular matrix holding all *squared inter-distances*
(zero on diagonal), using the
specified `metric`, of type [Metric::Enumerated type](@ref),
giving rise to distance function ``δ``. See [`distanceSqr`](@ref).

Only the lower triangular part is computed in order to optimize memory use.

By default, the result matrix is of type `Float32`. The type can be changed
to another real `type` using method (2).

<optional keyword arguments>:
- if ⏩=true (default) the computation of inter-distances is multi-threaded.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

**See**: [distance](@ref).

**See also**: [`laplacian`](@ref), [`laplacianEigenMaps`](@ref), [`spectralEmbedding`](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of 8 random 10x10 SPD matrices
Pset=randP(10, 8) # or, using unicode: 𝐏=randP(10, 8)
# Compute the squared inter-distance matrix according to the log Euclidean metric.
# This is much faster as compared to the Fisher metric and in general
# it is a good approximation.
Δ²=distanceSqrMat(logEuclidean, Pset)

# return a matrix of type Float64
Δ²64=distanceSqrMat(Float64, logEuclidean, Pset)

# Get the full matrix of inter-distances
fullΔ²=Hermitian(Δ², :L)
```
"""
function distanceSqrMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;
                         ⏩=true) where T<:AbstractFloat
   k, n, thr = dim(𝐏, 1), dim(𝐏, 2), nthreads()
   △=𝕃{type}(diagm(0 => zeros(k)))
   ⏩ && thr > 1 ? threaded=true : threaded=false
   #if threaded R=[(i, j) for j=1:k, i=1:k if i<j]; m=length(R) end # ranges
   if threaded
      R=[j for j=1:k, i=1:k if i<j]
      C=[i for j=1:k, i=1:k if i<j]
      m=length(R)
   end # ranges

   if metric == invEuclidean
       if threaded
           𝐏𝓲=ℍVector(undef, k)
           @threads for j=1:k 𝐏𝓲[j]=inv(𝐏[j]) end
           @threads for i=1:m △[R[i], C[i]]=ss(ℍ(𝐏𝓲[R[i]] - 𝐏𝓲[C[i]])) end
       else
           𝐏𝓲=map(inv, 𝐏) # inv preserves Hermitianity
           for j=1:k-1, i=j+1:k △[i, j]=ss(ℍ(𝐏𝓲[i] - 𝐏𝓲[j]))  end
       end

   elseif metric == logEuclidean
       if threaded
           𝐏𝓵=ℍVector(undef, k)
           @threads for j=1:k 𝐏𝓵[j]=ℍ(log(𝐏[j])) end
           @threads for i=1:m △[R[i], C[i]]=ss(ℍ(𝐏𝓵[R[i]] - 𝐏𝓵[C[i]])) end
       else
           𝐏𝓵=map(log, 𝐏)
           for j=1:k-1, i=j+1:k △[i, j]=ss(ℍ(𝐏𝓵[i] - 𝐏𝓵[j]))  end
       end

   elseif metric == ChoEuclidean
       if threaded
           𝐏L=𝕃Vector(undef, k)
           @threads for j=1:k 𝐏L[j]=choL(𝐏[j]) end
           @threads for i=1:m △[R[i], C[i]]=ss(𝐏L[R[i]] - 𝐏L[C[i]]) end
       else
           𝐏L=map(choL, 𝐏)
           for j=1:k-1, i=j+1:k △[i, j]=ss(𝐏L[i] - 𝐏L[j])  end
       end

   elseif metric==logCholesky
       if threaded
           𝐏L=𝕃Vector(undef, k)
           @threads for j=1:k 𝐏L[j]=choL(𝐏[j]) end
           @threads for i=1:m △[R[i], C[i]]=sst(tril(𝐏L[R[i]], -1)-tril(𝐏L[C[i]], -1), -1) + ssd(𝑓𝔻(log, 𝐏L[R[i]])-𝑓𝔻(log, 𝐏L[C[i]])) end
       else
           𝐏L=map(choL, 𝐏)
           for j=1:k-1, i=j+1:k △[i, j]=sst(tril(𝐏L[i], -1)-tril(𝐏L[j], -1), -1) + ssd(𝑓𝔻(log, 𝐏L[i])-𝑓𝔻(log, 𝐏L[j])) end
       end

   elseif metric==Jeffrey
       if threaded
           𝐏𝓲=ℍVector(undef, k)
           @threads for j=1:k 𝐏𝓲[j]=inv(𝐏[j]) end
           @threads for i=1:m △[R[i], C[i]]=0.5(tr(𝐏𝓲[C[i]], 𝐏[R[i]]) + tr(𝐏𝓲[R[i]], 𝐏[C[i]])) - n end
       else
           𝐏𝓲=map(inv, 𝐏)
           for j=1:k-1, i=j+1:k △[i, j]=0.5(tr(𝐏𝓲[j], 𝐏[i]) + tr(𝐏𝓲[i], 𝐏[j])) - n end
       end

   elseif metric==VonNeumann  # using formula: tr( PlogP + QLoqQ - PlogQ - QlogP)/2
       if threaded
           𝐏𝓵=ℍVector(undef, k)
           v=Vector(undef, k)
           @threads for j=1:k 𝐏𝓵[j]=ℍ(log(𝐏[j])); v[j]=tr(𝐏[j], 𝐏𝓵[j]) end
           @threads for i=1:m △[R[i], C[i]]=0.5*real(v[R[i]]+v[C[i]]-tr(𝐏[R[i]], 𝐏𝓵[C[i]])-tr(𝐏[C[i]], 𝐏𝓵[R[i]])) end
       else
           𝐏𝓵=[ℍ(log(P)) for P in 𝐏]
           v=[tr(𝐏[i], 𝐏𝓵[i]) for i=1:length(𝐏)]
           for j=1:k-1, i=j+1:k △[i, j]=0.5*real(v[i]+v[j]-tr(𝐏[i], 𝐏𝓵[j])-tr(𝐏[j], 𝐏𝓵[i])) end
       end

   elseif metric==Wasserstein
       if threaded
           𝐏½=ℍVector(undef, k)
           @threads for j=1:k 𝐏½[j]=sqrt(𝐏[j]) end
           @threads for i=1:m △[R[i], C[i]]=tr(𝐏[R[i]]) + tr(𝐏[C[i]]) -2tr(sqrt(ℍ(𝐏½[R[i]] * 𝐏[C[i]] * 𝐏½[R[i]]))) end
       else
           𝐏½=map(x->ℍ(sqrt(x)), 𝐏) # sqrt does not preserved hermitianity
           for j=1:k-1, i=j+1:k △[i, j]=tr(𝐏[i]) + tr(𝐏[j]) -2tr(sqrt(ℍ(𝐏½[i] * 𝐏[j] * 𝐏½[i]))) end
       end

   elseif metric ∈ (Euclidean, Fisher, logdet0)
       if threaded
           @threads for i=1:m △[R[i], C[i]]=distanceSqr(metric, 𝐏[R[i]], 𝐏[C[i]]) end
       else
           for j in 1:k-1, i in j+1:k △[i, j]=distanceSqr(metric, 𝐏[i], 𝐏[j])  end
       end

   else   @error("in RiemannianGeometryP.distanceSqrMat or .distanceMat function
                     (PosDefManifold Package): the chosen 'metric' does not exist")
   end # If metric

   return △
end #function

distanceSqrMat(metric::Metric, 𝐏::ℍVector;
                ⏩=true) = distanceSqrMat(Float32, metric, 𝐏; ⏩=⏩)

distance²Mat=distanceSqrMat


"""
```julia
    (1) distanceMat(metric::Metric, 𝐏::ℍVector;
    <⏩=true>)

    (2) distanceMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;
    <⏩=true>) where T<:AbstractFloat
```

Given a 1d array ``𝐏`` of ``k`` positive definite matrices
``{P_1,...,P_k}`` of [ℍVector type](@ref), create the ``k⋅k`` real
`LowerTriangular` matrix comprising elements
``δ(P_i, P_j)\\textrm{, for all }i>=j``.

This is the lower triangular matrix holding all *inter-distances*
(zero on diagonal), using the
specified `metric`, of type [Metric::Enumerated type](@ref),
giving rise to distance ``δ``. See [`distance`](@ref).

Only the lower triangular part is computed in order to optimize memory use.

By default, the result matrix is of type `Float32`. The type can be changed
to another real `type` using method (2).

The elements of this matrix are the square root of
[`distanceSqrMat`](@ref).

<optional keyword arguments>:
- if ⏩=true the computation of inter-distances is multi-threaded.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

**See**: [distance](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of 4 random 10x10 SPD matrices
Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
Δ=distanceMat(Fisher, Pset)

# return a matrix of type Float64
Δ64=distanceMat(Float64, Fisher, Pset)

# Get the full matrix of inter-distances
fullΔ=Hermitian(Δ, :L)
```
"""
distanceMat(type::Type{T}, metric::Metric, 𝐏::ℍVector;
            ⏩=true) where T<:AbstractFloat =
    sqrt.(distanceSqrMat(type, metric, 𝐏; ⏩=⏩))

distanceMat(metric::Metric, 𝐏::ℍVector;
            ⏩=true) = 
    sqrt.(distanceSqrMat(metric, 𝐏, ⏩=⏩))



"""
```julia
    laplacian(Δ²::𝕃{S}, epsilon::Real=0;
              <densityInvariant=false>) where S<:Real
```

Given a `LowerTriangular` matrix of squared inter-distances ``Δ^2``,
return the lower triangular part of the so-called
*normalized Laplacian* or *density-invariant normalized Laplacian*,
which in both cases is a symmetric Laplacian.
The elements of the Laplacian are of the same type as the elements of ``Δ^2``.
The result is a `LowerTriangular` matrix.

The definition of Laplacian given by Lafon (2004)[🎓](@ref) is implemented:

First, a [Gaussian radial basis functions](https://bit.ly/1HVyf55),
known as *Gaussian kernel* or *heat kernel*,
is applied to all elements of ``Δ^2``, such as

``W_{ij} = exp\\bigg(\\frac{\\displaystyle{-Δ^2_{ij}}}{\\displaystyle{2ε}}\\bigg)``,

where ``ε`` is the *bandwidth* of the kernel.

If <optional keyword argument> `densityInvariant=true` is used,
then the density-invariant transformation is applied

``W \\leftarrow E^{-1}WE^{-1}``

where ``E`` is the diagonal matrix holding on the main diagonal
the sum of the rows (or columns) of ``W``.

Finally, the normalized Laplacian (density-invariant or not) is defined as

``Ω = D^{-1/2}WD^{-1/2}``,

where ``D`` is the diagonal matrix holding on the main diagonal
the sum of the rows (or columns) of ``W``.

If you do not provide argument `epsilon`, the bandwidth ``ε`` is set to the
median of the elements of squared distance matrix ``Δ^2_{ij}``.
Another educated guess is the dimension of the original data, that is,
the data that has been used to compute the squared distance matrix.
For positive definite matrices this is ``n(n-1)/2``, where ``n`` is the
dimension of the matrices. Still another is the dimension of the ensuing
[`spectralEmbedding`](@ref) space.
Keep in mind that by tuning the `epsilon` parameter
(which must be positive) you can control both the rate of compression of the
embedding space and the spread of points in the embedding space.
See Coifman et *al.* (2008)[🎓](@ref) for a discussion on ``ε``.

!!! note "Nota Bene"
    The Laplacian as here defined can be requested for any
    input matrix of squared inter-distances, for example,
    those obtained on scalars or on vectors using appropriate metrics.
    In any case, only the lower triangular part of the Laplacian is
    taken as input. See [typecasting matrices](@ref).

**See also**: [`distanceSqrMat`](@ref), [`laplacianEigenMaps`](@ref), [`spectralEmbedding`](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of 4 random 10x10 SPD matrices
Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
Δ²=distanceSqrMat(Fisher, Pset)
Ω=laplacian(Δ²)

# density-invariant Laplacian
Ω=laplacian(Δ²; densityInvariant=true)

# increase the bandwidth
r=size(Δ², 1)
myεFactor=0.1
med=Statistics.median([Δ²[i, j] for j=1:r-1 for i=j+1:r])
ε=2*myεFactor*med
Ω=laplacian(Δ², ε; densityInvariant=true)
```
"""
 function laplacian(Δ²::𝕃{T}, epsilon::Real=0;
                    densityInvariant=false) where T<:Real

     # vector of sum of columns of W, using only the lower triangle Ω of W
     function sumCol(Ω::𝕃{T}, r::Int) where T<:Real
         D=zeros(T, r)
         for i=1:r
             for j=1:i   D[i] += Ω[i, j] end
             for l=i+1:r D[i] += Ω[l, i] end # conj(L[l, i]) for complex matrices
         end
         return D
     end

     # product DWD with D diagonal, using only the lower triangle Ω of W
     # overwrite Ω with the lower triangle of the result
     DWD!(D::Vector{T}, Ω::𝕃{T}, r::Int) where T<:Real =
                 for j=1:r, i=j:r Ω[i, j]*=D[i]*D[j] end

     r=size(Δ², 1)

     # Apply a Gaussian (heat) kernel
     epsilon≈0 ? ε=2*Statistics.median([Δ²[i, j] for j=1:r-1 for i=j+1:r]) : ε=2*epsilon
     Ω=𝕃{T}(diagm(0 => ones(r)))
     for j=1:r-1, i=j+1:r Ω[i, j]=exp(-Δ²[i, j]/ε) end

     if densityInvariant DWD!(inv.(sumCol(Ω, r)), Ω, r) end
     DWD!(sqrt.(inv.(sumCol(Ω, r))), Ω, r)
     return Ω # see laplacianEigenMaps
 end


"""
```julia
    laplacianEigenMaps(Ω::𝕃{S}, q::Int;
    <
    tol::Real=0.,
    maxiter::Int=300,
    verbose=false >) where S<:Real
```

**alias**: `laplacianEM`

Given the lower triangular part of a Laplacian ``Ω``
(see [`laplacian`](@ref) ) return the *eigen maps* in ``q`` dimensions,
i.e., the ``q`` eigenvectors of the Laplacian associated with the largest ``q``
eigenvalues, excluding the first (which is always equal to 1.0).
The eigenvectors are of the same type as ``Ω``. They are all divided
element-wise by the first eigenvector (see Lafon, 2004[🎓](@ref)).

The eigenvectors of the Laplacian are computed by the
power iterations+modified Gram-Schmidt method (see [`powerIterations`](@ref)),
allowing the execution of this function for Laplacian matrices of very large size.

Return the 4-tuple ``(Λ, U, iterations, convergence)``, where:
- ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps,
- ``U`` holds in columns the ``q`` eigenvectors holding the ``q`` coordinates of the points in the embedding space,
- ``iterations`` is the number of iterations executed by the power method,
- ``convergence`` is the convergence attained by the power method.

Using the notion of Laplacian, spectral embedding seek a
low-dimension representation of the data emphasizing local neighbothood
information while neglecting long-distance information.
The embedding is non-linear, however the embedding space is Euclidean.
The eigenvectors of ``U`` holds the coordinates of the points in the
embedding space (typically two- or three-dimensional for plotting or more
for clustering). Spectral embedding is done for plotting data in low-dimension,
clustering, imaging, classification, following their trajectories over time
or other dimensions, and much more.
For examples of applications see Ridrigues et *al.* (2018) [🎓](@ref)
and references therein.

**Arguments**:
- ``Ω`` is a real `LowerTriangular` normalized Laplacian obtained by the [`laplacian`](@ref) function,
- ``q`` is the dimension of the Laplacian eigen maps;
- The following are *<optional keyword arguments>* for the power iterations:
* `tol` is the tolerance for convergence (see below),
* `maxiter` is the maximum number of iterations allowed,
* if `verbose` is true, the convergence at all iterations will be printed.

!!! note "Nota Bene"
    The maximum value of ``q`` that can be requested is ``n-1``,
    where ``n`` is the size of the Laplacian.
    In general, ``q=2`` or ``q=3`` is requested.

    ``tol`` defaults to the square root of `Base.eps` of the (real) type
    of ``Ω``. This corresponds to requiring equality for the convergence criterion
    over two successive power iterations of about half of the significant digits.

**See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`spectralEmbedding`](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of 4 random 10x10 SPD matrices
Pset=randP(10, 4)
Δ²=distanceSqrMat(Fisher, Pset)
Ω=laplacian(Δ²)
evalues, maps, iterations, convergence=laplacianEM(Ω, 2)
evalues, maps, iterations, convergence=laplacianEM(Ω, 2; verbose=true)
evalues, maps, iterations, convergence=laplacianEM(Ω, 2; verbose=true, maxiter=500)
```

"""
function laplacianEigenMaps(Ω::𝕃{T}, q::Int;
                            tol::Real=0.,
                            maxiter::Int=300,
                            verbose=false) where T<:Real
    # make a check for q<size(Ω, 1)
    tol≈0. ? tolerance = √eps(T) : tolerance = tol
    (Λ, U, iter, conv) = powIter(Ω, q+1;
                                evalues=true,
                                tol=tolerance,
                                maxiter=maxiter,
                                verbose=verbose)
    for i=2:q+1 U[:, i]=U[:, i]./U[:, 1] end
    return 𝔻(Λ[2:q+1, 2:q+1]), U[1:size(U, 1), 2:q+1], iter, conv
end
laplacianEM=laplacianEigenMaps


"""
```julia
    (1) spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int, epsilon::Real=0;
    <
    tol::Real=0.,
    maxiter::Int=300,
    densityInvariant=false,
    verbose=false,
    ⏩=true >)

    (2) spectralEmbedding(type::Type{T}, metric::Metric, 𝐏::ℍVector, q::Int, epsilon::Real=0;
    < same optional keyword arguments as in (1) >) where T<:Real
```

**alias**: `spEmb`

Given a 1d array ``𝐏`` of ``k`` positive definite matrices ``{P_1,...,P_k}``
(real or complex), compute its *eigen maps* in ``q`` dimensions.

This function runs one after the other the functions:
- [`distanceSqrMat`](@ref) (compute the squared inter-distance matrix),
- [`laplacian`](@ref) (compute the normalized Laplacian),
- [`laplacianEigenMaps`](@ref) (get the eigen maps).

By default all computations above are done with `Float32` precision.
Another real type can be requested using method (2), where the `type` argument
is defined.

Return the 4-tuple `(Λ, U, iterations, convergence)`, where:
- ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps,
- ``U`` holds in columns the ``q`` eigenvectors holding the ``q`` coordinates of the points in the embedding space,
- ``iterations`` is the number of iterations executed by the power method,
- ``convergence`` is the convergence attained by the power method.

**Arguments**:
- `metric` is the metric of type [Metric::Enumerated type](@ref) used for computing the inter-distances,
- ``𝐏`` is a 1d array of ``k`` positive matrices of [ℍVector type](@ref),
- ``q`` is the dimension of the Laplacian eigen maps,
- ``epsilon`` is the bandwidth of the Laplacian (see [`laplacian`](@ref));
- The following *<optional keyword argument>* applyies for computing the inter-distances:
* if `⏩=true` (default) the computation of inter-distances is multi-threaded.
- The following *<optional keyword argument>* applyies to the computation of the Laplacian by the [`laplacian`](@ref) function:
* if `densityInvariant=true` the density-invariant Laplacian is computed (see [`laplacian`](@ref)).
- The following are *<optional keyword arguments>* for the power method iterative algorithm invoked by [`laplacianEigenMaps`](@ref):
    * `tol` is the tolerance for convergence of the power method (see below),
    * `maxiter` is the maximum number of iterations allowed for the power method,
    * if `verbose=true` the convergence at all iterations will be printed;

!!! note "Nota Bene"
    ``tol`` defaults to the square root of `Base.eps` of the `Float32` type (1)
    or of the `type` passed as argumant (2). This corresponds to requiring
    equality for the convergence criterion over two successive power iterations
    of about half of the significant digits.

    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

**See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`laplacianEigenMaps`](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of k random 10x10 SPD matrices
k=10
Pset=randP(10, k)
evalues, maps, iter, conv=spectralEmbedding(Fisher, Pset, 2)

# show convergence information
evalues, maps, iter, conv=spectralEmbedding(Fisher, Pset, 2; verbose=true)

# use Float64 precision.
evalues, maps, iter, conv=spectralEmbedding(Float64, Fisher, Pset, 2)

using Plots
# check eigevalues and eigenvectors
plot(diag(evalues))
plot(maps[:, 1])
plot!(maps[:, 2])
plot!(maps[:, 3])

# plot the data in the embedded space
plot(maps[:, 1], maps[:, 2], seriestype=:scatter, title="Spectral Embedding", label="Pset")

# try a different value of epsilon
evalues, maps, iter, conv=spEmb(Fisher, Pset, k-1, 0.01; maxiter=1000)
plot(maps[:, 1], maps[:, 2], seriestype=:scatter, title="Spectral Embedding", label="Pset")
# see the example in `Laplacian` function for more on this
```
"""
function spectralEmbedding(type::Type{T}, metric::Metric, 𝐏::ℍVector, q::Int, epsilon::Real=0;
                           tol::Real=0.,
                           maxiter::Int=300,
                           densityInvariant=false,
                           verbose=false,
                           ⏩=true) where T<:Real

    tol≈0. ? tolerance = √eps(type) : tolerance = tol
    return (Λ, U, iter, conv) =
            laplacianEM(laplacian(distance²Mat(type, metric, 𝐏, ⏩=⏩), epsilon;
                        densityInvariant=densityInvariant), q;
                        tol=tolerance,
                        maxiter=maxiter,
                        verbose=verbose)
end

spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int, epsilon::Real=0;
                  tol::Real=0.,
                  maxiter::Int=300,
                  densityInvariant=false,
                  verbose=false,
                  ⏩=true) =
    spectralEmbedding(Float32, metric, 𝐏, q, epsilon;
                      tol=tol,
                      maxiter=maxiter,
                      densityInvariant=densityInvariant,
                      verbose=verbose,
                      ⏩=⏩)

spEmb=spectralEmbedding


# -----------------------------------------------------------
# 4. Means (centers of mass, barycenters, ...)
# -----------------------------------------------------------

"""
```julia
    (1) mean(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex

    (2) mean(metric::Metric, D::𝔻{T}, E::𝔻{T}) where T<:Real

    (3) mean(metric::Metric, 𝐏::ℍVector;
        <
        w::Vector=[],
        ✓w=true,
        init::Union{ℍ, Nothing}=nothing,
        tol::Real=0.,
        verbose=false,
        ⏩=true >)

    (4) mean(metric::Metric, 𝐃::𝔻Vector;
    < same optional keyword arguments as in (3) >)
```

(1) Mean of two positive definite matrices, passed in arbitrary order as
arguments ``P`` and ``Q``, using the specified `metric` of type
[Metric::Enumerated type](@ref). The order is arbitrary as all metrics
implemented in **PosDefManifold** are symmetric.
This is the midpoint of the geodesic.
For the weighted mean of two positive definite matrices use instead
the [`geodesic`](@ref) function.
``P`` and ``Q`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

(2) Like in (1), but for two real diagonal positive definite matrices
``D`` and ``E``.

(3) [Fréchet mean](@ref) of an 1d array ``𝐏`` of ``k`` positive definite
matrices ``𝐏={P_1,...,P_k}`` of [ℍVector type](@ref),
with optional non-negative real weights ``w={w_1,...,w_k}`` and using the
specified `metric`as in (1).

(4) [Fréchet mean](@ref) of an 1d array ``𝐃`` of ``k`` positive definite
matrices ``𝐃={D_1,...,D_k}`` of [𝔻Vector type](@ref),
with optional non-negative real weights ``w={w_1,...,w_k}`` and using the
specified `metric`as in (1).

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

Adopting the `Fisher`, `logdet0` and `Wasserstein` metric in (3) and the
`logdet0` metric in (4), the mean is computed by means of an iterative
algorithm. A particular initialization for these algorithms can be
provided passing an Hermitian matrix as *<optional keyword argument>* `init`.
The convergence for these algorithm is required with a tolerance
given by *<optional keyword argument>* `tol`.
if `verbose=true` the covergence attained at each iteration is printed.
Other information such as if the algorithm has diverged is also printed.
For more options in computing these means call directly
functions [`geometricMean`](@ref), [`logdet0Mean`](@ref)
and [`wasMean`](@ref), which are called hereby.
For the meaning of the `tol` default value see the documentation of
these functions. See also the robust mean function [`geometricpMean`](@ref),
which cannot be called from here. Notice that arguments `init` and `tol`
have an effect only for the aferomentioned metrics in methods (3) and (4).

For (3) and (4), if `⏩=true` (default),
the computation of the mean is multi-threaded for all metrics.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

## Math

The Fréchet mean of a set of ``k`` matrices ``{P_1, P_2,..., P_k}`` weighted by
``{w_1, w_2,..., w_k}:\\sum_{i=1}^{k}w_i=1`` for the supported metrics are,
for those with closed form expression:

| Metric   | weighted Fréchet mean |
|:----------:|:----------- |
|Euclidean | ``\\sum_{i=1}^{k}w_i P_i`` |
|invEuclidean| ``\\big(\\sum_{i=1}^{k}w_i P_i^{-1}\\big)^{-1}``|
|ChoEuclidean| ``TT^*``, where ``T=bL_P + aL_Q`` |
|logEuclidean| ``\\textrm{exp}\\big(\\sum_{i=1}^{k}w_i\\hspace{1pt} \\textrm{log}P_i \\big)``|
|logCholesky| ``TT^*``, where `` T=\\sum_{i=1}^{k}(w_kS_k)+\\sum_{i=1}^{k}(w_k\\textrm{log}D_k)``|
|Jeffrey | ``A^{1/2}\\big(A^{-1/2}HA^{-1/2}\\big)^{1/2}A^{1/2}`` |

 and for those that are found by an iterative algorithm and that verify an equation:

| Metric   | equation verified by the weighted Fréchet mean |
|:----------:|:----------- |
|Fisher | ``\\sum_{i=1}^{k}w_i\\textrm{log}\\big(G^{-1/2} P_k G^{-1/2}\\big)=0.``|
|logdet0| ``\\sum_{i=1}^{k}w_i\\big(\\frac{1}{2}P_i+\\frac{1}{2}G\\big)^{-1}=G^{-1}`` |
|VonNeumann | N.A.|
|Wasserstein| ``G=\\sum_{i=1}^{k}w_i\\big( G^{1/2}  P_i G^{1/2}\\big)^{1/2}`` |

**legend:** ``L_X``, ``S_X`` and ``D_X``
are the Cholesky lower triangle of ``X``, its strictly lower triangular part
and diagonal part, respectively (hence, ``S_X+D_X=L_X``,  ``L_XL_X^*=X``).
``A`` and ``H`` are the weighted arithmetic and weighted harmonic mean, respectively.

**See**: [geodesic](@ref), [mean](@ref), [Fréchet mean](@ref).

**Examples**
```julia
using LinearAlgebra, Statistics, PosDefManifold
# Generate 2 random 3x3 SPD matrices
P=randP(3)
Q=randP(3)
M=mean(logdet0, P, Q) # (1)
M=mean(Euclidean, P, Q) # (1)

# passing several matrices and associated weights listing them
# weights vector, does not need to be normalized
R=randP(3)
mean(Fisher, ℍVector([P, Q, R]); w=[1, 2, 3])

# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4)
weights=[1, 2, 3, 1]
# passing a vector of Hermitian matrices (ℍVector type)
M=mean(Euclidean, Pset; w=weights) # (2) weighted Euclidean mean
M=mean(Wasserstein, Pset)  # (2) unweighted Wassertein mean
# display convergence information when using an iterative algorithm
M=mean(Fisher, Pset; verbose=true)

# run multi-threaded when the number of matrices is high
using BenchmarkTools
Pset=randP(20, 160)
@benchmark(mean(logEuclidean, Pset; ⏩=false)) # single-threaded
@benchmark(mean(logEuclidean, Pset)) # multi-threaded
```
"""
mean(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex =
    if metric==Fisher && size(P, 1)>=120 #(faster proc: Congedo et al., 2005)
        λ, B=eigen(P, Q) # the eigenvalues of Q are all 1.0 after diagonalization
        A=(Diagonal(λ))^(0.25)*inv(B)
        return ℍ(A'*A)
    else
        return geodesic(metric, P, Q, 0.5)
    end

mean(metric::Metric, D::𝔻{T}, E::𝔻{T}) where T<:Real = geodesic(metric, D, E, 0.5)

function mean(metric::Metric, 𝐏::ℍVector;
              w::Vector=[],
              ✓w=true,
              init::Union{ℍ, Nothing}=nothing,
              tol::Real=0.,
              verbose=false,
              ⏩=true)

    # iterative solutions
    if  metric == Fisher
        (G, iter, conv) = gMean(𝐏; w=w, ✓w=✓w, init=init, tol=tol, verbose=verbose, ⏩=⏩);
        return G
    end

    if  metric == logdet0
        (G, iter, conv) = ld0Mean(𝐏; w=w, ✓w=✓w, init=init, tol=tol, verbose=verbose, ⏩=⏩);
        return G
    end

    if  metric == Wasserstein
        (G, iter, conv) = wasMean(𝐏; w=w, ✓w=✓w, init=init, tol=tol, verbose=verbose, ⏩=⏩);
        return G
    end

    # closed-form expressions and exit
    k, n, thr = dim(𝐏, 1), dim(𝐏, 2), nthreads()
    ⏩ && thr > 1 ? threaded=true : threaded=false
    threaded && metric == logCholesky ? 𝐐 = 𝕃Vector(undef, k) : nothing
    isempty(w) ? v=[] : v = _getWeights(w, ✓w)

    if  metric == Euclidean
        if threaded
            isempty(w) ? (return fVec(𝛍, 𝐏)) : (return fVec(𝚺, 𝐏; w=v))
        else
            isempty(w) ? (return ℍ(𝛍(𝐏))) : (return ℍ(𝚺(map(*, v, 𝐏))))
        end

    elseif metric == invEuclidean
        if threaded
            if isempty(w) return inv(fVec(𝛍, inv, 𝐏))
            else          return inv(fVec(𝚺, inv, 𝐏; w=v)) end
        else
            if isempty(w) return inv(ℍ(𝛍(inv, 𝐏)))
            else          return inv(ℍ(𝚺(map(*, v, map(inv, 𝐏))))) end
        end

    elseif metric == logEuclidean
        if threaded
            if isempty(w) return ℍ(exp(fVec(𝛍, log, 𝐏)))
            else          return ℍ(exp(fVec(𝚺, log, 𝐏; w=v))) end
        else
            if isempty(w) return ℍ(exp(ℍ(𝛍(log, 𝐏))))
            else          return ℍ(exp(ℍ(𝚺(map(*, v, map(log, 𝐏)))))) end
        end

    elseif metric == ChoEuclidean
        if threaded
            if isempty(w) L=fVec(𝛍, choL, 𝐏)
            else          L=fVec(𝚺, choL, 𝐏; w=v) end
        else
            isempty(w) ? L = 𝛍(choL, 𝐏) : L = 𝚺(map(*, v, map(choL, 𝐏)))
        end
        return ℍ(L*L')

    elseif metric == logCholesky
        if threaded      @threads for i=1:k 𝐐[i] = choL(𝐏[i]) end
        else             𝐐=map(choL, 𝐏) end

        if isempty(w)
            Z=𝛍(tril(L,-1) for L in 𝐐) + exp(𝛍(𝑓𝔻(log, L) for L in 𝐐))
        else
            Z=𝚺(ω*tril(L,-1) for (ω, L) in zip(v, 𝐐)) + exp(𝚺(ω*𝑓𝔻(log, L) for (ω, L) in zip(v, 𝐐)))
        end
        return ℍ(Z*Z')

    elseif metric == Jeffrey # can be further optimized. See Faraki et al., 2015
        return mean(Fisher, mean(Euclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩), mean(invEuclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩))

    elseif metric == VonNeumann
        @warn "function RiemannianGeometryP.mean and .geodesic not defined for metric $metric"

    else
        @error "in RiemannianGeometryP.mean function: the chosen 'metric' does not exist"
    end # if metric
end # function


function mean(metric::Metric, 𝐃::𝔻Vector;
              w::Vector=[],
              ✓w=true,
              init::Union{ℍ, Nothing}=nothing,
              tol::Real=0.,
              verbose=false,
              ⏩=true)

    # iterative solutions
    if metric == logdet0
        (G, iter, conv) = ld0Mean(𝐃; w=w, ✓w=✓w, init=init, tol=tol, verbose=verbose, ⏩=⏩); return G
    end

    # closed-form expressions and exit
    k, n, thr = dim(𝐃, 1), dim(𝐃, 2), nthreads()
    ⏩ && thr > 1 ? threaded=true : threaded=false
    isempty(w) ? v=[] : v = _getWeights(w, ✓w)

    if     metric == Euclidean
        if threaded
            if isempty(w) return fVec(𝛍, 𝐃) else return fVec(𝚺, 𝐃; w=v) end
        else
            if isempty(w) return 𝛍(𝐃) else return 𝚺(map(*, v, 𝐃)) end
        end

    elseif metric == invEuclidean
        if threaded
            if isempty(w) return inv(fVec(𝛍, inv, 𝐃))
            else          return inv(fVec(𝚺, inv, 𝐃; w=v)) end
        else
            if isempty(w) return inv(𝛍(inv, 𝐃))
            else          return inv(𝚺(map(*, v, map(inv, 𝐃)))) end
        end

    elseif metric in (logEuclidean, Fisher)
        if threaded
            if isempty(w) return exp(fVec(𝛍, log, 𝐃))
            else          return exp(fVec(𝚺, log, 𝐃; w=v)) end
        else
            if isempty(w) return exp(𝛍(log, 𝐃))
            else          return exp(𝚺(map(*, v, map(log, 𝐃)))) end
        end

    elseif metric == ChoEuclidean
        if threaded
            if isempty(w) L=fVec(𝛍, √, 𝐃)
            else          L=fVec(𝚺, √, 𝐃; w=v) end
        else
            isempty(w) ? L = 𝛍(√, 𝐃) : L = 𝚺(map(*, v, map(√, 𝐃)))
        end
        return L*L

    elseif metric == logCholesky
        if threaded
            if isempty(w) return exp((fVec(𝛍, log, 𝔻Vector(map(√, 𝐃)))))^2
            else          return exp((fVec(𝚺, log, 𝔻Vector(map(√, 𝐃)); w=v)))^2 end
        else
            if isempty(w) return exp((𝛍(log, map(√, 𝐃))))^2
            else          return exp((𝚺(map(*, v, map(log, map(√, 𝐃))))))^2 end
        end

    elseif metric == Jeffrey
        D=mean(Euclidean, 𝐃; w=w, ✓w=✓w, ⏩=⏩)
        return D*((inv(D)*mean(invEuclidean, 𝐃; w=w, ✓w=✓w, ⏩=⏩))^0.5)

    elseif metric == VonNeumann
        @warn "function RiemannianGeometry.mean and .geodesic not defined for metric $metric"

    elseif  metric == Wasserstein
        return generalizedMean(𝐃, 0.5; w=w, ✓w=✓w, ⏩=⏩)

    else
        @error "in RiemannianGeometry.mean function: the chosen 'metric' does not exist"
    end # if metric
end # function


"""
```julia
    (1) means(metric::Metric, 𝒫::ℍVector₂;
    <⏩=true>)

    (2) means(metric::Metric, 𝒟::𝔻Vector₂;
    <⏩=true>)
```

(1) Given a 2d array ``𝒫`` of positive definite matrices as an [ℍVector₂ type](@ref)
compute the [Fréchet mean](@ref) for as many [ℍVector type](@ref) objects
as hold in ``𝒫``, using the specified `metric` of type
[Metric::Enumerated type](@ref).
Return the means in a vector of `Hermitian` matrices, that is, as an `ℍVector` type.

(2) Given a 2d array ``𝒟`` of real positive definite matrices as an [𝔻Vector₂ type](@ref)
compute the [Fréchet mean](@ref) for as many [𝔻Vector type](@ref) objects
as hold in ``𝒟``, using the specified `metric` of type
[Metric::Enumerated type](@ref).
Return the means in a vector of `Diagonal` matrices, that is, as a `𝔻Vector` type.

The weigted Fréchet mean is not supported in this function.

If `⏩=true` (default) the computation of the means
is multi-threaded.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

**See also**: [`mean`](@ref).

**Examples**
```julia
using PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)
# Generate a set of 40 random 4x4 SPD matrices
Qset=randP(3, 40) # or, using unicode: 𝐐=randP(3, 40)
# listing directly ℍVector objects
means(logEuclidean, ℍVector₂([Pset, Qset])) # or: means(logEuclidean, ℍVector₂([𝐏, 𝐐]))
# note that [𝐏, 𝐐] is actually a ℍVector₂ type object

# creating and passing an object of ℍVector₂ type
sets=ℍVector₂(undef, 2) # or: 𝒫=ℍVector₂(undef, 2)
sets[1]=Pset # or: 𝒫[1]=𝐏
sets[2]=Qset # or: 𝒫[2]=𝐐
means(logEuclidean, sets) # or: means(logEuclidean, 𝒫)

# going multi-threated

# first, create 20 sets of 200 50x50 SPD matrices
sets=ℍVector₂([randP(50, 200) for i=1:20])

# How much computing time we save ?
# (example min time obtained with 4 threads & 4 BLAS threads)
using BenchmarkTools

# non multi-threaded, mean with closed-form solution
@benchmark(means(logEuclidean, sets; ⏩=false)) # (6.196 s)

# multi-threaded, mean with closed-form solution
@benchmark(means(logEuclidean, sets)) # (1.897 s)

sets=ℍVector₂([randP(10, 200) for i=1:10])

# non multi-threaded, mean with iterative solution
# wait a bit
@benchmark(means(Fisher, sets; ⏩=false)) # (4.672 s )

# multi-threaded, mean with iterative solution
@benchmark(means(Fisher, sets)) # (1.510 s)
```

"""
means(metric::Metric, 𝒫::ℍVector₂; ⏩=true) =
        ℍVector([mean(metric, 𝐏; ⏩=⏩) for 𝐏 in 𝒫])

means(metric::Metric, 𝒟::𝔻Vector₂; ⏩=true) =
        𝔻Vector([mean(metric, 𝐃; ⏩=⏩) for 𝐃 in 𝒟])



"""
```julia
    generalizedMean(𝐏::Union{ℍVector, 𝔻Vector}, p::Real;
    <
    w::Vector=[],
    ✓w=true,
    ⏩=true >)
```

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices of
[ℍVector type](@ref) or real positive definite diagonal matrices of
[𝔻Vector type](@ref) and optional non-negative real weights vector
``w={w_1,...,w_k}``, return the *weighted generalized means* ``G``
with real parameter ``p``, that is,

``G=\\big(\\sum_{i=1}^{k}w_iP_i^p\\big)^{1/p}``.

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted generalized mean*

``G=\\big(\\sum_{i=1}^{k}P_i^p\\big)^{1/p}``.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.
This option is provided to allow
calling this function repeatedly without normalizing the weights each time.

If *<optional key argmuent>* ⏩=true the computation of the generalized mean
is multi-threaded.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

The following special cases for parameter ``p`` are noteworthy:
- For ``p=\\frac{1}{2}`` the generalized mean is the [modified Bhattacharyya mean](@ref).
- For ``p=1`` the generalized mean is the [Euclidean](@ref) mean.
- For ``p=-1`` the generalized mean is the [inverse Euclidean](@ref) mean.
- For (the limit of) ``p=0`` the generalized mean is the [log Euclidean](@ref) mean, which is the [Fisher](@ref) mean when matrices in 𝐏 all pair-wise commute.

Notice that when matrices in 𝐏 all pair-wise commute, for instance if the
matrices are diagonal,
the generalized means coincide with the [power means](@ref)
for any ``p∈[-1, 1]`` and for ``p=0.5`` it coincides also with the
[Wasserstein](@ref) mean. For this reason the generalized means are used
as default initialization of both the [`powerMean`](@ref) and [`wasMean`](@ref)
algorithm.

**See**: [generalized means](@ref).

**See also**: [`powerMean`](@ref), [`wasMean`](@ref), [`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, Statistics, PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)

# weights vector, does not need to be normalized
weights=[1, 2, 3, 1]

# unweighted mean
G = generalizedMean(Pset, 0.25) # or: G = generalizedMean(𝐏, 0.25)

# weighted mean
G = generalizedMean(Pset, 0.5; w=weights)

# with weights previously normalized we can set ✓w=false
weights=weights./sum(weights)
G = generalizedMean(Pset, 0.5; w=weights, ✓w=false)

# run multi-threaded when the number of matrices is high
using BenchmarkTools
Pset=randP(20, 160)
@benchmark(generalizedMean(Pset; ⏩=false)) # single-threaded
@benchmark(generalizedMean(Pset)) # multi-threaded
```
"""
function generalizedMean(𝐏::Union{ℍVector, 𝔻Vector}, p::Real;
                         w::Vector=[],
                         ✓w=true,
                         ⏩=true)
    𝕋=typeofMatrix(𝐏)
    if     p ≈ -1 return mean(invEuclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩)
    elseif p ≈  0 return mean(logEuclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩)
    elseif p ≈  1 return mean(Euclidean, 𝐏;    w=w, ✓w=✓w, ⏩=⏩)
    else
        k, n, thr = dim(𝐏, 1), dim(𝐏, 2), nthreads()
        ⏩ && thr > 1 ? threaded=true : threaded=false
        isempty(w) ? v=[] : v = _getWeights(w, ✓w)

        if threaded
            if isempty(w) return (fVec(𝛍, x->x^p, 𝐏))^(1/p)
            else          return (fVec(𝚺, x->x^p, 𝐏; w=v))^(1/p) end
        else
            if isempty(w) return 𝕋(𝛍(P^p for P in 𝐏))^(1/p)
            else          return 𝕋(𝚺(ω*P^p for (ω, P) in zip(v, 𝐏)))^(1/p) end
        end
    end # if p
end # function


"""
```julia
    geometricMean(𝐏::Union{ℍVector, 𝔻Vector};
    <
    w::Vector=[],
    ✓w=true,
    init=nothing,
    tol::Real=0.,
    maxiter::Int=200,
    adaptStepSize::Bool=true,
    verbose=false,
    ⏩=true >)
```

**alias**: `gMean`

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices of
[ℍVector type](@ref) or diagonal matrices of [𝔻Vector type](@ref)
and optional non-negative real weights vector ``w={w_1,...,w_k}``,
return the 3-tuple ``(G, iter, conv)``, where ``G`` is the mean according
to the [Fisher](@ref) metric and ``iter``, ``conv`` are the number of iterations
and convergence attained by the algorithm.
Mean ``G`` is the unique positive definite matrix satisfying

``\\sum_{i=1}^{k}w_i\\textrm{log}\\big(G^{-1/2} P_i G^{-1/2}\\big)=0.``

For estimating it, this function implements the well-known gradient descent
algorithm, but with an exponential decaying step size ``ς``, yielding iterations

``G ← G^{1/2}\\textrm{exp}\\big(ς\\sum_{i=1}^{k}w_i\\textrm{log}(G^{-1/2} P_i G^{-1/2})\\big)G^{1/2}.``

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted geometric mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

The following are more *<optional keyword arguments*>:
- `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the [Euclidean](@ref) mean will be used,
- `tol` is the tolerance for the convergence (see below).
- `maxiter` is the maximum number of iterations allowed
- if `verbose`=true, the convergence attained at each iteration and the step size ``ς`` is printed. Also, a *warning* is printed if convergence is not attained.
- if ⏩=true the iterations are multi-threaded (see below).
- if `adaptStepSize`=false the step size `ς` is fixed to 1 at all iterations.

If the input is a 1d array of ``k`` real positive definite diagonal matrices
the solution is available in closed-form as the log Euclidean
mean, hence the *<optional keyword arguments*> `init`, `tol` and `verbose`
have no effect and return the 3-tuple ``(G, 1, 0)``.
See the [log Euclidean](@ref) metric.

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges and `verbose` is true a **warning** is printed
    indicating the iteration when this happened.

    The exponential decaying step size features a faster convergence rate
    as compared to the fixed step size ``ς=1`` that is usually adopted.
    The decaying rate is inversely proportional to `maxiter`, thus,
    increase/decrease `maxiter` in order to set a slower/faster
    decaying rate. `maxiter` should not be set too low though.

    ``tol`` defaults to the square root of `Base.eps` of the nearest
    real type of data input ``𝐏``. This corresponds to requiring the
    norm of the satisfying matrix equation divided by the number of elements
    to vanish for about half the significant digits.

**See**: [Fisher](@ref) metric.

**See also**: [`geometricpMean`](@ref), [`powerMean`](@ref),
 [`wasMean`](@ref), [`logdet0Mean`](@ref), [`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)

# unweighted mean
G, iter, conv = geometricMean(Pset) # or G, iter, conv = geometricMean(𝐏)

# weights vector, does not need to be normalized
weights=[1, 2, 3, 1]

# weighted mean
G, iter, conv = geometricMean(Pset, w=weights)

# print the convergence at all iterations
G, iter, conv = geometricMean(Pset; verbose=true)

# now suppose Pset has changed a bit, initialize with G to hasten convergence
Pset[1]=ℍ(Pset[1]+(randP(3)/100))
G, iter, conv = geometricMean(Pset; w=weights, ✓w=true, verbose=true, init=G)

# run multi-threaded when the number of matrices is high
using BenchmarkTools
Pset=randP(20, 120)
@benchmark(geometricMean(Pset; ⏩=false)) # single-threaded
@benchmark(geometricMean(Pset)) # multi-threaded

# show the mean and the input points using spectral embedding
using Plots
k=80
Pset=randP(20, k)
G, iter, conv = geometricMean(Pset)
push!(Pset, G)
Λ, U, iter, conv=spectralEmbedding(Fisher, Pset, 2; verbose=true)
plot(U[1:k, 1], U[1:k, 2], seriestype=:scatter, title="Spectral Embedding", label="Pset")
plot!(U[k+1:k+1, 1], U[k+1:k+1, 2], seriestype=:scatter, label="mean")
```
"""
function geometricMean( 𝐏::ℍVector;
                        w::Vector=[],
                        ✓w=true,
                        init=nothing,
                        tol::Real=0.,
                        maxiter::Int=200,
                        adaptStepSize::Bool=true,
                        verbose=false,
                        ⏩=true)

    (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v) = _setVar_IterAlg(𝐏, w, ✓w, tol, ⏩)
    _giveStartInfo_IterAlg(threaded, verbose, "geometricMean Fixed-Point")
    init === nothing ? M = mean(Euclidean, 𝐏; w=v, ✓w=false, ⏩=⏩) : M = ℍ(init)
    💡 = similar(M, type) # new iteration solution
    if threaded 𝐐 = 𝕄Vector(repeat([𝐏[1]], thr)) end # memory pre-allocation for fVec function
    c1(M⁻½::ℍ, 𝐏::ℍVector) = cong(M⁻½, 𝐏, ℍVector) # utility function
    c2(M⁻½::ℍ, P::ℍ) = cong(M⁻½, P, ℍ) # utility function

    # M -< M½ { exp[ς( w_i{sum(i=1 to k) log(M⁻½ 𝐏[i] M⁻½)} )] } M½
    while true
        M½, M⁻½ = pow(M, 0.5, -0.5)
        if threaded
            isempty(w) ? ∇ = fVec(𝛍, log, c1(M⁻½, 𝐏), allocs=𝐐) :
                         ∇ = fVec(𝚺, log, c1(M⁻½, 𝐏), w=v, ✓w=false, allocs=𝐐)
        else
            isempty(w) ? ∇ = ℍ(𝛍(log(c2(M⁻½, P)) for P in 𝐏)) :
                         ∇ = ℍ(𝚺(ω * log(c2(M⁻½, P)) for (ω, P) in zip(v, 𝐏)))
        end
        adaptStepSize ? 💡 = ℍ(M½*exp(ς*∇)*M½) : 💡 = ℍ(M½*exp(∇)*M½)

        conv = norm(∇)/n² # norm of the satisfying equation. It must vanish upon convergence
        if adaptStepSize ς = exp(-ℯ * golden * iter / maxiter) end # exponetially decaying step size
        verbose && println("iteration: ", iter, "; convergence: ", conv, "; ς: ", round(ς * 1000)/1000)
        (diverging = conv > oldconv) && verbose && @warn("geometricMean diverged at:", iter)
        (overRun = iter == maxiter) && @warn("geometricMean reached the max number of iterations before convergence:", iter)
        (converged = conv <= tolerance) || overRun==true ? break : M = 💡
        oldconv=conv; iter += 1
    end # while

    _giveEndInfo_IterAlg(converged, verbose)
    return (💡, iter, conv)
end



geometricMean(𝐃::𝔻Vector;
              w::Vector=[],
              ✓w=true,
              init=nothing,
              tol::Real=0.,
              verbose=false,
              ⏩=true) = mean(logEuclidean, 𝐃; w=w, ✓w=false, ⏩=⏩), 1, 0

gMean=geometricMean



"""
```julia
    geometricpMean(𝐏::ℍVector, p::Real=0.5;
    <
    w::Vector=[],
    ✓w=true,
    init=nothing,
    tol::Real=0.,
    maxiter::Int=500,
    adaptStepSize=true,
    verbose=false,
    ⏩=true >)
```

**alias**: `gpmean`

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices of
[ℍVector type](@ref), a real parameter ``0<p<=1`` and optional non-negative real
weights vector ``w={w_1,...,w_k}``, return the 3-tuple ``(G, iter, conv)``,
where ``G`` is the *p-mean*, i.e., the mean according to the
[Fisher](@ref) metric minimizing the *p-dispersion* (see below) and
``iter``, ``conv`` are the number of
iterations and convergence attained by the algorithm.

This function implements the p-dispersion gradient descent
algorithm with step-size ``ς`` (to be published), yielding iterations

``G ←G^{1/2}\\textrm{exp}\\big(ς\\sum_{i=1}^{k}pδ^2(G, P_i)^{p-1}w_i\\textrm{log}(G^{-1/2} P_i G^{-1/2})\\big)G^{1/2}.``

- if ``p=1`` this yields the geometric mean (implemented specifically in [`geometricMean`](@ref)).
- if ``p=0.5`` this yields the geometric median (default).

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted geometric-p mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

The following are more *<optional keyword arguments*>:
* `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the [Euclidean](@ref) mean will be used,
* `tol` is the tolerance for the convergence (see below).
* `maxiter` is the maximum number of iterations allowed.
* if `adaptStepSize`=true (default) the step size ``ς`` for the gradient descent is adapted at each iteration (see below).
* if `verbose`=true, the step-size and convergence attained at each iteration is printed. Also, a *warning* is printed if convergence is not attained.
* if ⏩=true the iterations are multi-threaded (see below).

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

    If the algorithm diverges and `verbose` is true a **warning** is printed
    indicating the iteration when this happened. This algorithm may temporary
    diverge, still reach convergence. Overall, while all other iterative
    algorithms implemented in **PosDefMaifold** are very stable, this is not.

    The smaller the parameter ``p`` is, the slower and less likely the
    convergence is. If the algorithm does not converge, try increasing ``p``,
    initializing the algorithm with the output of [`geometricMean`](@ref)
    and/or eliminating the otliers from the input set ``𝐏``.

    If `adaptStepSize` is true (default) the step-size ``ς`` is adapted at
    each iteration, otherwise a fixed step size ``ς=1`` is used.
    Adapting the step size in general hastens convergence and improves
    the convergence behavior.

    ``tol`` defaults to the square root of `Base.eps` of the nearest
    real type of data input ``𝐏``. This corresponds to requiring the
    norm of the satisfying matrix equation divided by the number of elements
    to vanish for about half the significant digits.

**See**: [Fisher](@ref) metric.

**See also**: [`geometricMean`](@ref), [`powerMean`](@ref),
[`wasMean`](@ref), [`logdet0Mean`](@ref), [`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold, Plots

# This examples show that this algorithm is more robust to outliers
# as compared to the standard geometric mean algorithm

# Generate a set of 100 random 10x10 SPD matrices
Pset=randP(10, 100)

# Get the usual geometric mean for comparison
G, iter1, conv1 = geometricMean(Pset, verbose=true)

# change p to observe how the convergence behavior changes accordingly
# Get the median (default)
H, iter2, conv2 = geometricpMean(Pset, verbose=true)
# Get the p-mean for p=0.25
H, iter2, conv2 = geometricpMean(Pset, 0.25, verbose=true)

println(iter1, " ", iter2); println(conv1, " ", conv2)

# move the first matrix in Pset to possibly create an otlier
Pset[1]=geodesic(Fisher, G, Pset[1], 3)
G1, iter1, conv1 = geometricMean(Pset, verbose=true)
H1, iter2, conv2 = geometricpMean(Pset, 0.25, verbose=true)
println(iter1, " ", iter2); println(conv1, " ", conv2)

# collect the geometric and p-means, before and after the
# introduction of the outier in vector of Hermitian matrices `S`.
S=HermitianVector([G, G1, H, H1])

# check the interdistance matrix Δ² to verify that the geometric mean
# after the introduction of the outlier (``G1``) is farther away from
# the geometric mean as compared to how much ``H1`` is further away
# from ``H``, i.e., that element (4,3) is much smaller than element (2,1).
Δ²=distance²Mat(Float64, Fisher, S)

# how far are all these matrices from all the others?
fullΔ²=Hermitian(Δ², :L)
dist=[sum(fullΔ²[:, i]) for i=1:size(fullΔ², 1)]

# plot the matrices in `S` using spectral embedding.
using Plots
Λ, U, iter, conv = laplacianEM(laplacian(Δ²), 3;  verbose=true)
plot([U[1, 1]], [U[1, 2]], seriestype=:scatter, label="g-mean")
plot!([U[2, 1]], [U[2, 2]], seriestype=:scatter, label="g-mean outlier")
plot!([U[3, 1]], [U[3, 2]], seriestype=:scatter, label="p-mean")
plot!([U[4, 1]], [U[4, 2]], seriestype=:scatter, label="p-mean outlier")

# estimate how much you gain running the algorithm in multi-threaded mode
using BenchmarkTools
Pset=randP(20, 120)
@benchmark(geometricpMean(Pset; ⏩=true)) # single-threaded
@benchmark(geometricpMean(Pset)) # multi-threaded
```
"""
function geometricpMean(𝐏::ℍVector, p::Real=goldeninv;
                        w::Vector = [], ✓w = true,
                        init = nothing,
                        tol::Real = 0.,
                        maxiter::Int = 500,
                        adaptStepSize = true,
                        verbose = false,
                        ⏩= false)

    (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v) = _setVar_IterAlg(𝐏, w, ✓w, tol, ⏩)
    _giveStartInfo_IterAlg(threaded, verbose, "geometricpMean Fixed-Point")
    𝑓, d², q, ςHasNotChanged, ςold = Fisher, distance², p-1, 0, 0
    init === nothing ? M = mean(Euclidean, 𝐏; w=v, ✓w=false, ⏩=⏩) : M = ℍ(init)
    💡 = similar(M, type)
    𝐑 = similar(𝐏)
    if threaded 𝐐 = similar(𝐏) end

    while true
        M½, M⁻½ = pow(M, 0.5, -0.5)
        #M -< M^1/2 {  exp[epsilon( 1/n{sum(i=1 to n) ln(M^-1/2 Mi M^-1/2)} )] } M^1/2
        if threaded
            @threads for i=1:k 𝐑[i] = ℍ(M⁻½ * 𝐏[i] * M⁻½) end
            if isempty(w)
                @threads for i=1:k 𝐐[i] = p * d²(𝑓, 𝐑[i])^q * log(𝐑[i]) end
                ∇ = fVec(𝛍, 𝐐)
            else
                @threads for i=1:k 𝐐[i] =  p * d²(𝑓, 𝐑[i])^q * v[i] * log(𝐑[i]) end
                ∇ = fVec(𝚺, 𝐐)
            end
        else
            for i=1:k 𝐑[i] = ℍ(M⁻½ * 𝐏[i] * M⁻½) end
            if isempty(w)
                ∇ = ℍ(𝛍(p * d²(𝑓, R)^q * log(R) for R in 𝐑))
            else
                ∇ = ℍ(𝚺(p * d²(𝑓, R)^q * ω * log(R) for (ω, R) in zip(v, 𝐑)))
            end
        end

        conv = norm(∇) / n²

        if adaptStepSize
            if  conv<=oldconv
                💡 = ℍ(M½ * exp(ς*∇) * M½)
                ς < 1 ? ς/=1.1 : ς*=1.1
                oldconv=conv
            else
                ς=1
            end
            #oldconv=conv
            verbose && println("iteration: ", iter, "; convergence: ", conv, "; ς: ", round(ς*1000)/1000)
        else
            💡 = ℍ(M½ * exp(∇) * M½)
            oldconv=conv
            verbose && println("iteration: ", iter, "; convergence: ", conv)
        end

        (diverging = conv > oldconv) && verbose && @warn("geometricpMean diverged at:", iter)
        (overRun = iter == maxiter) && @warn("geometricpMean reached the max number of iterations before convergence:", iter)
        (converged = conv <= tolerance) || ς <= tolerance || overRun==true ? break : M = 💡
        iter += 1

    end # while

    _giveEndInfo_IterAlg(converged, verbose)
    return (💡, iter, conv)
end

gpMean=geometricpMean


"""
```julia
    logdet0Mean(𝐏::Union{ℍVector, 𝔻Vector};
    <
    w::Vector=[],
    ✓w=true,
    init=nothing,
    tol::Real=0.,
    maxiter::Int=500,
    verbose=false,
    ⏩=true >)
```

**alias**: `ld0Mean`

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices of
[ℍVector type](@ref) or real positive definite diagonal matrices of
[𝔻Vector type](@ref) and optional
non-negative real weights vector ``w={w_1,...,w_k}``,
return the 3-tuple ``(G, iter, conv)``, where ``G`` is the mean according
to the [logdet zero](@ref) metric and ``iter``, ``conv`` are the number of iterations
and convergence attained by the algorithm.
Mean ``G`` is the unique positive definite matrix satisfying

``\\sum_{i=1}^{k}w_i\\big(\\frac{1}{2}P_i+\\frac{1}{2}G\\big)^{-1}-G^{-1}=0``.

For estimating it, this function implements the fixed-point iteration algorithm
suggested by (Moakher, 2012, p315)[🎓](@ref), yielding iterations

``G ← \\frac{1}{2}\\big(\\sum_{i=1}^{k}w_i(P_i+G)^{-1}\\big)^{-1}``.

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted logdet zero mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

The following are more *<optional keyword arguments*>:
- `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the [Euclidean](@ref) mean will be used,
- `tol` is the tolerance for the convergence (see below).
- `maxiter` is the maximum number of iterations allowed.
- if `verbose`=true, the convergence attained at each iteration is printed and a *warning* is printed if convergence is not attained.
- if ⏩=true the iterations are multi-threaded (see below).

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges and `verbose` is true a **warning** is printed
    indicating the iteration when this happened.

    ``tol`` defaults to 100 times the square root of `Base.eps` of the nearest
    real type of data input ``𝐏``. This corresponds to requiring the
    square root of the relative
    convergence criterion over two successive iterations to vanish for about
    half the significant digits minus 2.

**See**: [logdet zero](@ref) metric, [modified Bhattacharyya mean](@ref).

**See also**: [`powerMean`](@ref), [`wasMean`](@ref), [`logdet0Mean`](@ref),
[`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)

# unweighted mean
G, iter, conv = logdet0Mean(Pset) # or G, iter, conv = logdet0Mean(𝐏)

# weights vector, does not need to be normalized
weights=[1, 2, 3, 1]

# weighted mean
G, iter, conv = logdet0Mean(Pset, w=weights)

# print the convergence at all iterations
G, iter, conv = logdet0Mean(Pset; w=weights, verbose=true)

# suppose Pset has changed a bit; initialize with G to hasten convergence
Pset[1]=ℍ(Pset[1]+(randP(3)/100))
G, iter, conv = logdet0Mean(Pset; w=weights, ✓w=false, verbose=true, init=G)

# estimate how much you gain running the algorithm in multi-threaded mode
using BenchmarkTools
Pset=randP(20, 120)
@benchmark(logdet0Mean(Pset; ⏩=false)) # single-threaded
@benchmark(logdet0Mean(Pset)) # multi-threaded
```
"""
function logdet0Mean(𝐏::Union{ℍVector, 𝔻Vector};
                    w::Vector=[],
                    ✓w=true,
                    init=nothing,
                    tol::Real=0.,
                    maxiter::Int = 500,
                    verbose=false,
                    ⏩=true)


    (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v) = _setVar_IterAlg(𝐏, w, ✓w, tol, ⏩)
    _giveStartInfo_IterAlg(threaded, verbose, "logDet0Mean Fixed-Point")
    𝕋, l = typeofMatrix(𝐏), k/2
    init === nothing ? M = mean(Euclidean, 𝐏; w=v, ✓w=false, ⏩=⏩) : M = 𝕋(init)
    💡 = similar(M, type)
    if threaded 𝐐 = similar(𝐏) end

    while true
        if threaded
            if isempty(w)
                @threads for i=1:k 𝐐[i] = inv(𝕋(𝐏[i]+M)) end
                💡 = l * inv(fVec(𝚺, 𝐐))
            else
                @threads for i=1:k 𝐐[i] = v[i] * inv(𝕋(𝐏[i]+M)) end
                💡 = 0.5 * inv(fVec(𝚺, 𝐐))
            end
        else
            if isempty(w)
                💡 = l * inv(𝕋(𝚺(inv(𝕋(P+M)) for P in 𝐏)))
            else
                💡 = 0.5 * inv(𝕋(𝚺(ω * inv(𝕋(P+M)) for (ω, P) in zip(v, 𝐏))))
            end
        end

        conv = norm(💡-M)/n²
        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (diverging = conv > oldconv) && verbose && @warn("logdet0Mean diverged at:", iter)
        (overRun = iter == maxiter) && @warn("logdet0Mean reached the max number of iterations before convergence:", iter)
        (converged = conv <= tolerance) || overRun==true ? break : M = 💡
        oldconv=conv
        iter += 1
    end # while

    _giveEndInfo_IterAlg(converged, verbose)
    return (💡, iter, conv)
end

ld0Mean=logdet0Mean


"""
```julia
    wasMean(𝐏::Union{ℍVector, 𝔻Vector};
    <
    w::Vector=[],
    ✓w=true,
    init=nothing,
    tol::Real=0.,
    maxiter::Int=500,
    verbose=false,
    ⏩=true >)
```

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices
of [ℍVector type](@ref) or real positive definite diagonal matrices of
[𝔻Vector type](@ref) and optional non-negative real weights vector
``w={w_1,...,w_k}``,
return the 3-tuple ``(G, iter, conv)``, where ``G`` is the mean according
to the [Wasserstein](@ref) metric and ``iter``, ``conv`` are the number of iterations
and convergence attained by the algorithm.
Mean ``G`` is the unique positive definite matrix satisfying

``G=\\sum_{i=1}^{k}w_i\\big( G^{1/2}  P_i G^{1/2}\\big)^{1/2}``.

For estimating it, this function implements the fixed-point iterative algorithm
proposed by (Álvarez-Esteban et *al.*, 2016)[🎓](@ref):

``G ← G^{-1/2}\\big(\\sum_{i=1}^{k} w_i(G^{1/2}P_i G^{1/2})^{1/2}\\big)^2 G^{-1/2}``.

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted Wassertein mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and they should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

The following are more *<optional keyword arguments*>:
- `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with ``p=0.5`` will be used,
- `tol` is the tolerance for the convergence (see below).
- `maxiter` is the maximum number of iterations allowed.
- if `verbose`=true, the convergence attained at each iteration is printed and a *warning* is printed if convergence is not attained.
- if ⏩=true the iterations are multi-threaded (see below).

If the input is a 1d array of ``k`` real positive definite diagonal matrices
the solution is available in closed-form as the modified Bhattacharyya mean,
hence the *<optional keyword arguments*> `init`, `tol` and `verbose`
have no effect and return the 3-tuple ``(G, 1, 0)``.
See [modified Bhattacharyya mean](@ref).

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges and `verbose` is true a **warning** is printed indicating
    the iteration when this happened.

    ``tol`` defaults to the square root of `Base.eps` of the nearest
    real type of data input ``𝐏``. This corresponds to requiring the
    norm of the satisfying matrix equation divided by the number of elements
    to vanish for about half the significant digits.

**See**: [Wasserstein](@ref) metric.

**See also**: [`powerMean`](@ref), [`wasMean`](@ref), [`logdet0Mean`](@ref),
[`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)

# unweighted mean
G, iter, conv = wasMean(Pset) # or: G, iter, conv = wasMean(𝐏)

# weights vector, does not need to be normalized
weights=[1, 2, 3, 1]

# weighted mean
G, iter, conv = wasMean(Pset; w=weights)

# print the convergence at all iterations
G, iter, conv = wasMean(Pset; w=weights, verbose=true)

# suppose 𝐏 has changed a bit; initialize with G to hasten convergence
Pset[1]=ℍ(Pset[1]+(randP(3)/100))
G, iter, conv = wasMean(Pset; w=weights, verbose=true, init=G)

# estimate how much you gain running the algorithm in multi-threaded mode
using BenchmarkTools
Pset=randP(20, 120)
@benchmark(wasMean(Pset; ⏩=false)) # single-threaded
@benchmark(wasMean(Pset)) # multi-threaded
```
"""
function wasMean(𝐏::ℍVector;
                w::Vector=[],
                ✓w=true,
                init=nothing,
                tol::Real=0.,
                maxiter::Int = 500,
                verbose=false,
                ⏩=true)

    (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v) = _setVar_IterAlg(𝐏, w, ✓w, tol, ⏩)
    _giveStartInfo_IterAlg(threaded, verbose, "wasMean Fixed-Point")
    init === nothing ? M = generalizedMean(𝐏, 0.5; w=v, ✓w=false, ⏩=⏩) : M = ℍ(init)
    💡 = similar(M, type)
    if threaded 𝐐 = similar(𝐏) end

    while true
        M½, M⁻½ = pow(M, 0.5, -0.5)
        if threaded
            if isempty(w)
                @threads for i=1:k 𝐐[i] = ℍ(√ℍ(M½*𝐏[i]*M½)) end
                💡 = ℍ(M⁻½ * sqr(fVec(𝛍, 𝐐)) * M⁻½)
            else
                @threads for i=1:k 𝐐[i] = v[i] * ℍ(√ℍ(M½*𝐏[i]*M½)) end
                💡 = ℍ(M⁻½ * sqr(fVec(𝚺, 𝐐)) * M⁻½)
            end
        else
            if isempty(w)
                💡 = ℍ(M⁻½ * sqr(ℍ(𝛍(√ℍ(M½*P*M½) for P in 𝐏))) * M⁻½)
            else
                💡 = ℍ(M⁻½ * sqr(ℍ(𝚺((√ℍ(M½*P*M½) * ω) for (ω, P) in zip(v, 𝐏)))) * M⁻½)
            end
        end

        #conv = √norm(💡-M)/norm(M)
        conv = norm(💡-M)/n²
        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (diverging = conv > oldconv) && verbose && @warn("wasMean diverged at:", iter)
        (overRun = iter == maxiter) && @warn("wasMean reached the max number of iterations before convergence:", iter)
        (converged = conv <= tolerance) || overRun==true ? break : M = 💡
        oldconv=conv
        iter += 1
    end # while

    _giveEndInfo_IterAlg(converged, verbose)
    return (💡, iter, conv)
end

wasMean(𝐃::𝔻Vector;
        w::Vector=[], ✓w=true, init=nothing, tol::Real=0., verbose=false, ⏩=true) =
        generalizedMean(𝐃, 0.5, w=w, ✓w=✓w, ⏩=⏩), 1, 0


"""
```julia
    powerMean(𝐏::Union{ℍVector, 𝔻Vector}, p::Real;
    <
    w::Vector=[],
    ✓w=true,
    init=nothing,
    tol::Real=0.,
    maxiter::Int=500,
    verbose=false,
    ⏩=true >)
```

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` positive definite matrices
of [ℍVector type](@ref) or real positive definite diagonal matrices of
[𝔻Vector type](@ref), an optional non-negative real weights vector
``w={w_1,...,w_k}`` and a real parameter `p` ``\\in[-1, 1]``, return the
3-tuple ``(G, iter, conv)``, where ``G`` is
Lim and Palfia (2012)'s [power means](@ref)  of order ``p`` and
``iter``, ``conv`` are the number of iterations
and convergence attained by the algorithm, respectively.
Mean ``G`` is the unique positive definite matrix satisfying

``G=\\sum_{i=1}^{k}w_i(G\\textrm{#}_pP_i)``,

where ``G\\textrm{#}_pP_i`` is the [Fisher](@ref) [geodesic](@ref) equation.
In particular:

- with ``p=-1`` this is the *harmonic mean* (see the [inverse Euclidean](@ref) metric),
- with ``p=+1`` this is the *arithmetic mean* (see the [Euclidean](@ref) metric),
- at the limit of ``p`` evaluated at zero from both side this is the *geometric mean* (see [Fisher](@ref) metric).

For estimating power means for ``p\\in(-1, 1)``, this function implements
the  fixed-point iterative algorithm of (Congedo et *al.*, 2017b)[🎓](@ref).
For ``p=0`` (geometric mean)
this algorithm is run two times with a small positive and negative value
of ``p`` and the geometric mean of the two
resulting means is returned, as suggested in (Congedo et *al.*, 2017b)[🎓](@ref).
This way of estimating the geometric mean of
a set of matrices is faster as compared to the usual gradient descent algorithm.

If you don't pass a weight vector with *<optional keyword argument>* ``w``,
return the *unweighted power mean*.

If *<optional keyword argument>* `✓w=true` (default), the weights are
normalized so as to sum up to 1, otherwise they are used as they are passed
and should be already normalized.  This option is provided to allow
calling this function repeatedly without normalizing the same weights
vector each time.

The following are more *<optional keyword arguments*>:
- `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with parameter ``p`` will be used.
- `tol` is the tolerance for the convergence (see below).
- `maxiter` is the maximum number of iterations allowed.
- if `verbose`=true, the convergence attained at each iteration is printed and a *warning* is printed if convergence is not attained.
- if ⏩=true the iterations are multi-threaded.

If the input is a 1d array of ``k`` real positive definite diagonal matrices
the solution is available in closed-form as the generalized
mean of order `p`, hence the *<optional keyword arguments*>
`init`, `tol` and `verbose`
have no effect and return the 3-tuple ``(G, 1, 0)``.
See [generalized means](@ref).

!!! note "Nota Bene"
    [Multi-threading](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1)
    is automatically disabled if Julia is instructed to use only one thread.
    See [Threads](@ref).

    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges and `verbose` is true a **warning** is printed indicating
    the iteration when this happened.

    ``tol`` defaults to the square root of `Base.eps` of the nearest
    real type of data input ``𝐏``. This corresponds to requiring the
    norm of the difference of the matrix solution over two successive
    iterations divided by the number of elements in the matrix
    to vanish for about half the significant digits.

(2) Like in (1), but for a 1d array ``𝐃={D_1,...,D_k}`` of ``k``
real positive definite diagonal matrices of [𝔻Vector type](@ref).
In this case the solution is available in closed-form, hence the
*<optional keyword arguments*> `init`, `tol` and `verbose` have no effect and return
the 3-tuple ``(G, 1, 0)``. See [generalized means](@ref).

**See**: [power means](@ref), [generalized means](@ref), [modified Bhattacharyya mean](@ref).

**See also**: [`generalizedMean`](@ref), [`wasMean`](@ref), [`logdet0Mean`](@ref),
[`mean`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
# Generate a set of 4 random 3x3 SPD matrices
Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)

# unweighted mean
G, iter, conv = powerMean(Pset, 0.5) # or G, iter, conv = powerMean(𝐏, 0.5)

# weights vector, does not need to be normalized
weights=[1, 2, 3, 1]

# weighted mean
G, iter, conv = powerMean(Pset, 0.5; w=weights)

# print the convergence at all iterations
G, iter, conv = powerMean(Pset, 0.5; w=weights, verbose=true)

# suppose 𝐏 has changed a bit; initialize with G to hasten convergence
Pset[1]=ℍ(Pset[1]+(randP(3)/100))
G, iter, conv = powerMean(Pset, 0.5; w=weights, verbose=true, init=G)

# estimate how much you gain running the algorithm in multi-threaded mode
using BenchmarkTools
Pset=randP(20, 120)
@benchmark(powerMean(Pset, 0.5; ⏩=false)) # single-threaded
@benchmark(powerMean(Pset, 0.5)) # multi-threaded
```
"""
function powerMean(𝐏::ℍVector, p::Real;
         w::Vector=[],
         ✓w=true,
         init=nothing,
         tol::Real=0.,
         maxiter::Int=500,
         verbose=false,
         ⏩=true)

  if ! (-1<=p<=1)
       @error("The parameter p for power means must be in range [-1...1]")
  else

    if p ≈-1
       return (mean(invEuclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩), 1, 0)
    elseif p ≈ 0
       LE=mean(logEuclidean, 𝐏, w=w, ✓w=✓w, ⏩=⏩)
       P, iter1, conv1=powerMean(𝐏,  0.01; w=w, ✓w=✓w, init=LE, tol=tol, verbose=verbose, ⏩=⏩)
       Q, iter2, conv2=powerMean(𝐏, -0.01; w=w, ✓w=✓w, init=P, tol=tol, verbose=verbose, ⏩=⏩)
       return (geodesic(Fisher, P, Q,  0.5), iter1+iter2, (conv1+conv2)/2)

    elseif p ≈ 1
       return (mean(Euclidean, 𝐏; w=w, ✓w=✓w, ⏩=⏩), 1, 0)
    else
       (k, n, type, thr, n², iter, conv, oldconv, converged, ς, threaded, tolerance, v) = _setVar_IterAlg(𝐏, w, ✓w, tol, ⏩)
       _giveStartInfo_IterAlg(threaded, verbose, "powerMean Fixed-Point")
       absp, sqrtn = abs(p), √n
       r = -0.375/absp
       init === nothing ? M = generalizedMean(𝐏, p; w=v, ✓w=false, ⏩=⏩) : M = ℍ(init)
       p<0 ? X=ℍ(M^(0.5)) : X=ℍ(M^(-0.5))
       💡, H, 𝒫 = similar(X, type), similar(X, type), similar(𝐏)
       p<0 ? 𝒫=[inv(P) for P in 𝐏] : 𝒫=𝐏
       if threaded 𝐐 = similar(𝐏) end

       while true
           if threaded
               if isempty(w)
                   @threads for i=1:k 𝐐[i] = ℍ(X*𝒫[i]*X')^absp end
                   H=fVec(𝛍, 𝐐)
               else
                   @threads for i=1:k 𝐐[i] = v[i] * ℍ(X*𝒫[i]*X')^absp end
                   H=fVec(𝚺, 𝐐)
               end
           else
               if isempty(w)
                   H=𝛍(ℍ(X*P*X')^absp for P in 𝒫)
               else
                   H=𝚺(ω * ℍ(X*P*X')^absp for (ω, P) in zip(v, 𝒫))
               end
           end
           💡 = ℍ(H)^r * X

       conv = norm(H-I)/n²
       verbose && println("iteration: ", iter, "; convergence: ", conv)
       (diverging = conv > oldconv) && verbose && @warn("powerMean diverged at:", iter)
       (overRun = iter == maxiter) && @warn("powerMean: reached the max number of iterations before convergence:", iter)
       (converged = conv <= tolerance) || overRun==true ? break : X = 💡
       oldconv=conv
       iter += 1
       end # while

    end # if

    _giveEndInfo_IterAlg(converged, verbose)
    p<0 ? (return ℍ((💡)'*💡), iter, conv) : (return inv(ℍ((💡)'*💡)), iter, conv)
  end # if !(-1<=p<=1)
end

powerMean(𝐃::𝔻Vector, p::Real;
          w::Vector=[],
          ✓w=true,
          init=nothing,
          tol::Real=0.,
          maxiter::Int=500,
          verbose=false,
          ⏩=true) = generalizedMean(𝐃, p; w=w, ✓w=✓w, ⏩=⏩), 1, 0



"""
```julia
    (1) inductiveMean(metric::Metric, 𝐏::ℍVector)

    (2) inductiveMean(metric::Metric, 𝐏::ℍVector, q::Int, Q::ℍ)
```

**alias**: `indMean`

(1) Compute the Fréchet mean of 1d array ``𝐏={P_1,...,P_k}`` of ``k``
positive definite matrices of [ℍVector type](@ref) with a law of large
number inductive procedure (Ho et *al.,* 2013; Lim and Palfia, 2019;
Massart et *al.*, 2018)[🎓](@ref), such as

``G_1=P_1,``

``G_i=γ(i^{-1}, G_{(i-1)}, P_i), i=2,...,k,``

where ``γ(i^{-1}, G_{(i-1)}, P_i)`` is a step on the [geodesic](@ref) relying
``G_{(i-1)}`` to ``P_i`` with arclength ``i^{-1}``
using the specified `metric`, of type [Metric::Enumerated type](@ref).

(2) Like (1), but for the set of matrices ``𝐐 ∪ 𝐏``,
where it is assumed knowledge only of the set ``𝐏``,
the mean of ``𝐐`` (Hermitian matrix argument `Q`) and the number of
matrices in ``𝐐`` (integer argument `q`).
This method can be used, for example, for updating a block on-line algorithm,
where ``𝐏`` is the incoming block, `Q` the previous mean estimation
and `q` the cumulative number of matrices on which the mean has been
computed on-line.

For Fréchet means that do not have a closed form expression,
this procedure features a computational complexity amounting to less than
two iterations of gradient descent or fixed-point algorithms. This comes at
the price of an approximation.
In fact, the solution is not invariant to permutations of the matrices
in array 𝐏 and convergence to the Fréchet mean for finite samples
is not ensured
(see Lim and Palfia, 2019; Massart et *al.*, 2018)[🎓](@ref).

Since the inductive mean uses the [`geodesic`](@ref) function,
it is not available for the Von Neumann metric.

**Examples**
```julia
# A set of 100 matrices for which we want to compute the mean
𝐏=randP(10, 100)

𝐏1=ℍVector(collect(𝐏[i] for i=1:50)) # first 50
𝐏2=ℍVector(collect(𝐏[i] for i=51:100)) # last 50

# inductive mean of the whole set 𝐏
G=inductiveMean(Fisher, 𝐏)

# mean using the usual gradient descent algorithm
H, iter, conv=geometricMean(𝐏)

# inductive mean of 𝐏 given only 𝐏2,
# the number of matrices in 𝐏1 and the mean of 𝐏1
G2=inductiveMean(Fisher, 𝐏2, length(𝐏1), mean(Fisher, 𝐏1))

# average error
norm(G-H)/(dim(G, 1)^2)
norm(G2-H)/(dim(G, 1)^2)
```
"""
function inductiveMean(metric::Metric, 𝐏::ℍVector)
    if metric ∉ (VonNeumann)
        G = 𝐏[1]
        for k=2:length(𝐏) G=geodesic(metric, G, 𝐏[k], 1/k) end
        return G
    else
        @warn("The inductive mean is not available for the Von Neumann metric")
    end
end


"""
```julia
    midrange(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
```

Midrange (average of extremal values) of positive definite matrices
``P`` and ``Q``. Only the Fisher metric is supported, allowing the so-called
*geometric midrange*. This has been defined in Mostajeran et *al.* (2019)
[🎓](@ref) as

``P * Q = \\frac{1}{\\sqrt{\\lambda_(min)}+\\sqrt{\\lambda_(max)}}\\Big(Q+\\sqrt{\\lambda_(min)*\\lambda_(max)}P\\Big)``,

where ``\\lambda_(min)`` and ``\\lambda_(max)`` are the extremal generalized
eigenvalues of ``P`` and ``Q``.

**Examples**
```julia
P=randP(3)
Q=randP(3)
M=midrange(Fisher, P, Q)
```

"""
function midrange(metric::Metric, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    if metric == Fisher
        λ=eigvals(𝕄(A), 𝕄(B))
        λmin=λ[1]
        λmax=λ[end]
        return (1/(√λmin+√λmax)) * (Q + P*√(λmin*λmax))
    else @warn "The matrix midrange is available only for the Fisher metric."
    end
end

function inductiveMean(metric::Metric, 𝐏::ℍVector, q::Int, Q::ℍ)
    if metric ∉ (VonNeumann)
        G = Q
        for k=q+1:q+length(𝐏) G=geodesic(metric, G, 𝐏[k-q], 1/k) end
        return G
    else
        @warn("The inductive mean is not available for the Von Neumann metric")
    end
end

indMean=inductiveMean


# -----------------------------------------------------------
# 5. Tangent Space Tools
# -----------------------------------------------------------

"""
```julia
    (1) logMap(metric::Metric, P::ℍ{T}, G::ℍ{T})

    (2) logMap(metric::Metric, 𝐏::ℍVector, G::ℍ{T})
    for all the above: where T<:RealOrComplex
```

(1) *Logaritmic Map:* map a positive definite matrix ``P`` from the SPD or
Hermitian manifold into the tangent space at base-point ``G``
using the [Fisher](@ref) metric.

``P`` and ``G`` must be flagged as `Hermitian`.
See [typecasting matrices](@ref).

The map is defined as

`` Log_G(P)=S=G^{1/2}\\textrm{log}\\big(G^{-1/2}PG^{-1/2}\\big)G^{1/2}``.

`metric` is a metric of type [Metric::Enumerated type](@ref).

The result is an `Hermitian` matrix.

(2) *Logarithmic map* (1) at base-point ``G`` at once for ``k`` positive definite
matrices in 1d array ``𝐏={P_1,...,P_k}`` of [ℍVector type](@ref).

The result is an `ℍVector`.

!!! note "Nota Bene"
    Currently only the [Fisher](@ref) metric is supported for
    tangent space operations.

The inverse operation is [`expMap`](@ref).

**See also**: [`vecP`](@ref), [`parallelTransport`](@ref).

**Examples**
```julia
using PosDefManifold
(1)
P=randP(3)
Q=randP(3)
metric=Fisher
G=mean(metric, P, Q)
# projecting P at the base point given by the geometric mean of P and Q
S=logMap(metric, P, G)

(2)
Pset=randP(3, 4)
# projecting all matrices in Pset at the base point given by their geometric mean.
Sset=logMap(Fisher, Pset, mean(Fisher, Pset))
```
"""
function logMap(metric::Metric, P::ℍ{T}, G::ℍ{T}) where T<:RealOrComplex
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * log(ℍ(G⁻½ * P * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.logMap function:
                 only the Fisher metric is supported for the logarithmic map."
    end
end

function logMap(metric::Metric, 𝐏::ℍVector, G::ℍ{T}) where T<:RealOrComplex
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍVector([ℍ(G½ * log(ℍ(G⁻½ * P * G⁻½)) * G½) for P in 𝐏])
    else @warn "in RiemannianGeometryP.logMap function:
                 only the Fisher metric is supported for the logarithmic map."
    end
end


"""
```julia
    (1) expMap(metric::Metric, S::ℍ{T}, G::ℍ{T})

    (2) expMap(metric::Metric, 𝐒::ℍVector, G::ℍ{T})
    for all the above: where T<:RealOrComplex
```

(1) *Exponential Map:* map a tangent vector (a matrix) ``S`` from the tangent
space at base-point ``G`` into the SPD or Hermitian manifold
(using the [Fisher](@ref) metric).

``S`` and ``G`` must be flagged as `Hermitian`.
See [typecasting matrices](@ref).

The map is defined as

``Exp_G(S)=P=G^{1/2}\\textrm{exp}\\big(G^{-1/2}SG^{-1/2}\\big)G^{1/2}``.

`metric` is a metric of type [Metric::Enumerated type](@ref).

The result is an `Hermitian` matrix.

(2) *Exponential map* (1) at base-point ``G`` at once for ``k`` tangent vectors
(matrices) in 1d array ``𝐒={S_1,...,S_k}`` of [ℍVector type](@ref).

The result is an `ℍVector`.

!!! note "Nota Bene"
    Currently only the [Fisher](@ref) metric is supported for tangent space operations.

The inverse operation is [`logMap`](@ref).

**Examples**
```julia
(1)
using PosDefManifold, LinearAlgebra
P=randP(3)
Q=randP(3)
G=mean(Fisher, P, Q)
# projecting P on the tangent space at the Fisher mean base point G
S=logMap(Fisher, P, G)
# projecting back onto the manifold
P2=expMap(Fisher, S, G)

(2)
Pset=randP(3, 4)
# projecting all matrices in Pset at the base point given by their geometric mean.
G=mean(Fisher, Pset)
Sset=logMap(Fisher, Pset, G)
# projecting back onto the manifold
Pset2=expMap(Fisher, Sset, G)
```
"""
function expMap(metric::Metric, S::ℍ{T}, G::ℍ{T}) where T<:RealOrComplex
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * exp(ℍ(G⁻½ * S * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.expMap function:
              only the Fisher metric is supported for the exponential map"
    end
end

function expMap(metric::Metric, 𝐒::ℍVector, G::ℍ{T}) where T<:RealOrComplex
#    k, thr = dim(𝐒, 1), nthreads()
#    ⏩ && thr > 1 ? threaded=true : threaded=false
#    if threaded 𝐐 = similar(𝐒)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍVector([ℍ(G½ * exp(ℍ(G⁻½ * S * G⁻½)) * G½) for S in 𝐒])
    else @warn "in RiemannianGeometryP.expMap function:
              only the Fisher metric is supported for the exponential map"
    end
end


"""
```julia
    vecP(S::Union{ℍ{T}, Symmetric{R}};
         range::UnitRange=1:size(S, 2)) where T<:RealOrComplex where R<:Real =
```

*Vectorize* a tangent vector (which is an `Hermitian` or `Symmetric` matrix)
``S``:  mat ↦ vec.

It gives weight ``1`` to diagonal elements and ``√2`` to off-diagonal elements
so as to preserve the norm (Barachant et *al.*, 201E)[🎓](@ref), such as

``∥S∥_F=∥vecP(S)∥_F``.

The result is a vector holding ``n(n+1)/2`` elements, where ``n``
is the size of ``S``.

``S`` must be flagged as `Hermitian` or `Symmetric`.
See [typecasting matrices](@ref).

The reverse operation is provided by [`matP`](@ref),
which always return an `Hermitian` matrix.

If an optional keyword argument `range` is provided,
the vectorization concerns only the rows (or columns,
since the input matrix is symmetric or Hermitian)
in the range. Note that in this case the operation
cannot be reverted by the [`matP`](@ref), that is,
in this case the matrix is 'stuck' in the tangent space.

**Examples**
```julia
using PosDefManifold
P=randP(3)
Q=randP(3)
G=mean(Fisher, P, Q)
# projecting P at the base point given by the geometric mean of P and Q
S=logMap(Fisher, P, G)
# vectorize S
v=vecP(S)
# vectorize onlt the first two columns of S
v=vecP(S; range=1:2)
```
"""
vecP(S::Union{ℍ{T}, Symmetric{R}};
     range::UnitRange=1:size(S, 2)) where T<:RealOrComplex where R<:Real =
     [(if i==j return S[i, j] else return (S[i, j])*sqrt2 end) for j=range for i=j:size(S, 1)]


"""
```julia
    matP(ς::Vector{T}) where T<:RealOrComplex
```

*Matrizize* a tangent vector (vector) ς :  vec -> mat.

This is the function reversing the [`vecP`](@ref) function,
thus the weighting applied therein is reversed as well.

If ``ς=vecP(S)`` and ``S`` is a ``n⋅n`` Hermitian or Symmetric matrix,
``ς``  is a tangent vector of size ``n(n+1)/2``.
The result of calling `matP(ς)` is then ``n⋅n`` matrix ``S``.
``S`` is always returned flagged as `Hermitian`.


**To Do**: This function may be rewritten more efficiently.

**Examples**
```julia
using PosDefManifold
P=randP(3)
Q=randP(3)
G=mean(Fishr, P, Q)
# projecting P at onto the tangent space at the Fisher mean base point
S=logMap(Fisher, P, G)
# vectorize S
v=vecP(S)
# Rotate the vector by an orthogonal matrix
n=Int(size(S, 1)*(size(S, 1)+1)/2)
U=randP(n)
z=U*v
# Get the point in the tangent space
S=matP(z)
```
"""
function matP(ς::Vector{T}) where T<:RealOrComplex
  n=Int((-1+√(1+8*length(ς)))/2) # Size of the matrix whose vectorization vector v has size length(v)
  S=Matrix{T}(undef, n, n)
  l=0
  for j in 1:n-1
    l=l+1
    @inbounds S[j, j]=ς[l]
    for i in j+1:n
      l=l+1
      @inbounds S[i, j]=invsqrt2*ς[l]
    end
  end
  @inbounds S[n, n]=ς[end]
  return ℍ(S, :L)
end

"""
```julia
    (1) parallelTransport(S::ℍ{T}, P::ℍ{T}, Q::ℍ{T})

    (2) parallelTransport(S::ℍ{T}, P::ℍ{T})

    (3) parallelTransport(𝐒::ℍVector, P::ℍ{T}, Q::ℍ{T})

    (4) parallelTransport(𝐒::ℍVector, P::ℍ{T})
    for all the above: where T<:RealOrComplex
```

**alias**: `pt`

(1) *Parallel transport* of tangent vector ``S`` (a matrix) lying on
the tangent space at base-point ``P`` to the tangent space at base-point ``Q``.

``S``, ``P`` and ``Q`` must all be `Hermitian` matrices.
Return an `Hermitian` matrix.
The transport is defined as:

``∥_{(P→Q)}(S)=\\big(QP^{-1}\\big)^{1/2}S\\big(QP^{-1}\\big)^{H/2}``.

If ``S`` is a positive definite matrix in the manifold (and not a tangent vector)
it will be 'trasported' from ``P`` to ``Q``, amounting to (Yair et *al.*, 2019[🎓](@ref))
- project ``S`` onto the tangent space at base-point ``P``,
- parallel transport it to the tangent space at base-point ``Q``,
- project it back onto the manifold at base-point ``Q``.

(2) *Parallel transport* as in (1), but to the tangent space at base-point
the *identity matrix*.

The transport reduces in this case to:

``∥_{(P→I)}(S)=P^{-1/2}SP^{-1/2}``.

(3) *Parallel transport* as in (1) at once for ``k`` tangent vectors
(matrices) in 1d array ``𝐒={S_1,...,S_k}`` of [ℍVector type](@ref).

(4) *Parallel transport* as in (2) at once for ``k`` tangent vectors
(matrices) in 1d array ``𝐒={S_1,...,S_k}`` of [ℍVector type](@ref).

!!! note "Nota Bene"
    Currently only the [Fisher](@ref) metric is supported for parallel transport.

**See also**: [`logMap`](@ref), [`expMap`](@ref), [`vecP`](@ref), [`matP`](@ref).

**Examples**
```julia
using PosDefManifold

(1)
P=randP(3)
Q=randP(3)
G=mean(Fisher, P, Q)

# i. projecting P onto the tangent space at base-point G
S=logMap(Fisher, P, G)
# ii. parallel transport S to the tangent space at base-point Q
S_=parallelTransport(S, G, Q)
# iii. projecting back into the manifold at base-point Q
P_=expMap(Fisher, S_, Q)

# i., ii. and iii. can be done simply by
PP_=parallelTransport(P, G, Q)
# check
P_≈PP_ ? println(" ⭐ ") : println(" ⛔ ")

(2)
P=randP(3)
Q=randP(3)
G=mean(Fisher, P, Q)
# transport to the tangent space at base-point the identity
PP_=parallelTransport(P, G)

(3)
Pset=randP(3, 4)
Q=randP(3)
G=mean(Fisher, Pset)
# trasport at once all matrices in Pset
Pset2=parallelTransport(Pset, G, Q)

(4)
Pset=randP(3, 4)
G=mean(Fisher, Pset)
# recenter all matrices so to have mean=I
Pset2=parallelTransport(Pset, G)
# check
mean(Fisher, Pset2) ≈ I ? println(" ⭐ ") : println(" ⛔ ")
```
"""
function parallelTransport(S::ℍ{T}, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    # this is the classical definition of parallel transport.
    # IT DOES NOT TRANSPORT MATRICES ALONG THE MANIFOLD
    # M = mean(Fisher, P, Q)
    # W = sqrt(M) * inv(P) # curious: if W=M*inv(P), return=M
    T<:Real ? W = real(√(Q * inv(P))) : W = √(Q * inv(P))
    return ℍ(W * S * W')
end

function parallelTransport(S::ℍ{T}, P::ℍ{T}) where T<:RealOrComplex
    W = invsqrt(P)
    return ℍ(W * S * W')
end

function parallelTransport(𝐒::ℍVector, P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    T<:Real ? W = real(√(Q * inv(P))) : W = √(Q * inv(P))
    return ℍVector([ℍ(W * S * W') for S in 𝐒])
end

function parallelTransport(𝐒::ℍVector, P::ℍ{T}) where T<:RealOrComplex
    W = invsqrt(P)
    return ℍVector([ℍ(W * S * W') for S in 𝐒])
end

pt=parallelTransport


# -----------------------------------------------------------
# 6. Procrustes Problems
# -----------------------------------------------------------

"""
```julia
    procrustes(P::ℍ{T}, Q::ℍ{T}, extremum="min") where T<:RealOrComplex
```

Given two positive definite matrices ``P`` and ``Q``,
return by default the solution of problem

``\\textrm{argmin}_Uδ(P,U^HQU)``,

where ``U`` varies over the set of unitary matrices and ``δ(.,.)`` is a
distance or divergence function.

``U^HQU`` is named in physics the *unitary orbit* of ``Q``.

If the argument `extremum` is passed as "max", it returns instead the solution of

``\\textrm{argmax}_Uδ(P,U^HQU)``.

``P`` and ``Q`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

As it has been shown in Bhatia and Congedo (2019)[🎓](@ref),
using each of the [Fisher](@ref), [logdet zero](@ref), [Wasserstein](@ref)
and the Kullback-Leibler divergence (see [logdet α](@ref)),
the best approximant to ``P`` from the unitary orbit of ``Q``
commutes with ``P`` and, surprisingly, has the same closed-form expression, namely

``U_Q^↓U_P^{↓H}`` for the argmin and ``U_Q^↑U_P^{↓H}`` for the argmax,

where ``U^↓`` denotes the eigenvector matrix of the subscript argument with
eigenvectors in columns sorted by *decreasing* order of corresponding eigenvalues and
``U^↑`` denotes the eigenvector matrix of the subscript argument with
eigenvectors in columns sorted by *increasing* order of corresponding eigenvalues.

The same solutions are known since a long time also by solving the extremal
problem here above using the [Euclidean](@ref) metric (Umeyama, 1988).

The generalized Procrustes problem

``\\textrm{argmin}_Usum_{i=1}^{k}δ(P_i,U^HQ_iU)``

can be solved using Julia package [Manopt](https://github.com/kellertuer/Manopt.jl).

**Examples**
```julia
using PosDefManifold
P=randP(3)
Q=randP(3)
# argmin problem
U=procrustes(P, Q)
# argmax problem
V=procrustes(P, Q, "max")
```
"""
function procrustes(P::ℍ{T}, Q::ℍ{T}, extremum="min") where T<:RealOrComplex
    Pup=eigvecs(P)
    Qup=eigvecs(Q)
    Pdown=reverse(Pup, dims=(2))
    if      extremum=="min"
            Qdown=reverse(Qup, dims=(2))
            return Qdown*Pdown'
    elseif  extremum=="max"
            return Qup*Pdown'
    else    @warn "in RiemannianGeometryP.procrustes: the argument 'extremum' is incorrect."
    end
end
