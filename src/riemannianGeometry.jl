#    Unit riemannianGeometry.jl, part of PosDefManifold Package for julia language
#    v 0.1.2 - last update 22th of April 2019
#
#    MIT License
#    Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
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
# 0. Internal Functions
#    By convention their name begin with underscore char
# -----------------------------------------------------------
_attributes(𝐏::ℍVector)=( size(𝐏[1], 1), length(𝐏))

# Given a non-negative weight vector normalize the weights so as to sum up to 1
# if ✓w == true and if they are not already normalized
function _getWeights(w::Vector, ✓w::Bool, k::Int)
    if ✓w==true
        s=𝚺(w)
        if s ≉  1.0 return w./s else return w end
    else return w
    end # if
end



# -----------------------------------------------------------
# 1. Geodesic Equations
# -----------------------------------------------------------

"""
    geodesic(metric::Metric, P::ℍ, Q::ℍ, a::Real)

 Move along the [geodesic](@ref) from point ``P`` to point ``Q``
 (two positive definite matrices) with *arclegth* ``0<=a<=1``,
 using the specified metric, of type [Metric::Enumerated type](@ref).

 For all metrics,
 - with ``a=0`` we stay at ``P``,
 - with ``a=1`` we move up to ``Q``,
 - with ``a=1/2`` we move to the mid-point of ``P`` and ``Q`` (mean).

 Using the Fisher metric argument `a` can be *any* real number, for instance:
 - with ``0<a<1`` we move toward ``Q`` (*attraction*),
 - with ``a>1`` we move over and beyond ``Q`` (*extrapolation*),
 - with ``a<0`` we move back away from Q (*repulsion*).

 Note that if ``Q=I``, the Fisher geodesic move is simply ``P^a``
 (no need to call this funtion then).

 For the [logdet zero](@ref) and [Jeffrey](@ref) metric no closed form expression
 for the geodesic is available (to the best of authors' knowledge),
 so in this case the geodesic is found as the weighted mean using [`mean(@ref)`].
 For the [Von Neumann](@ref) not even an expression for the mean is available,
 so in this case the geodesic is not provided and a *warning* is printed.

 ``P`` and ``Q`` must be flagged by julia as `Hermitian`.
 See [typecasting matrices](@ref).

 **Maths**

 For points ``P``, ``Q`` and arclength ``a``, letting ``b=1-a``,
 the geodesic equations for the supported metrics are:

| Metric   | geodesic equation |
|:----------:|:----------- |
|Euclidean| ``bP + aQ`` |
|invEuclidean| ``\\big(bP^{-1} + aQ^{-1}\\big)^{-1}``|
|ChoEuclidean| ``TT^*``, where ``T=bL_P + aL_Q``|
|logEuclidean| ``\\text{exp}\\big(b\\text{log}(P) + a\\text{log}(Q)\\big)``|
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

 ## Examples
    using PosDefManifold
    P=randP(10)
    Q=randP(10)
    # Wasserstein mean
    M=geodesic(Wasserstein, P, Q, 0.5)
    # extrapolate suing the Fisher metric
    E=geodesic(Fisher, P, Q, 2)

"""
function geodesic(metric::Metric, P::ℍ, Q::ℍ, a::Real)
    if a ≈ 0 return P end
    if a ≈ 1 return Q end
    b = 1-a

    if      metric==Euclidean
    return  P*b + Q*a

    elseif  metric==invEuclidean
    return  inv( ℍ(inv(P)*b + inv(Q)*a) )

    elseif  metric==logEuclidean
    return  ℍ( exp( ℍ(log(P)*b + log(Q)*a) ) )

    elseif  metric==Fisher
            P½, P⁻½ = pow(P, 0.5, -0.5)
    return  ℍ( P½ * (P⁻½ * Q * P⁻½)^a * P½ )

    elseif  metric in (logdet0, Jeffrey)
    return  mean(metric, ℍVector([P, Q]), w=[b, a], ✓w=false)

    elseif  metric==VonNeumann
            @warn("An expression for the geodesic is not available for the Von neumann metric")

    elseif  metric==ChoEuclidean
            T=choL(P)*b + choL(Q)*a
    return  ℍ(T*T')

    elseif  metric==logCholesky
            LP=choL(P)
            LQ=choL(Q)
            slLP=tril(LP,-1)
            T=slLP+a*(tril(LQ,-1)-slLP)+𝑓𝑫(x->x, LP)*exp(a*(𝑓𝑫(log, LQ)-𝑓𝑫(log, LP)))
    return  ℍ(T*T')

    elseif  metric==Wasserstein
            if isreal(P) && isreal(Q)
                    return ℍ( (b^2)*P + (a^2)*Q + (a*b)*real(√(P*Q)+√(Q*P)) )
            else    return ℍ( (b^2)*P + (a^2)*Q + (a*b)*(√(P*Q)+√(Q*P)) )
            end

    else    @warn("in RiemannianGeometryP.geodesic function
                 (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end # function



# -----------------------------------------------------------
# 2. Distances
# -----------------------------------------------------------

"""
    (1) distanceSqr(metric::Metric, P::ℍ)
    (2) distanceSqr(metric::Metric, P::ℍ, Q::ℍ)

 **alias**: `distance²`

 (1) Return ``δ^2(P, I)``, the *square of the distance* (or *divergence*) of positive definite
 matrix ``P`` from the the identity matrix. See [distance from the origin](@ref).

 (2) Return ``δ^2(P, Q)``, the *square of the distance* (or *divergence*) between two
 positive definite matrices ``P`` and ``Q``. See [distance](@ref).

 In both cases the distance function ``δ`` is induced by the argument `metric` of type
 [Metric::Enumerated type](@ref).

 ``P`` in (1) and ``P``, ``Q`` in (2) must be flagged by julia as `Hermitian`.
 See [typecasting matrices](@ref).

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

 ## Examples (1)
    using PosDefManifold
    P=randP(10)
    d=distanceSqr(Wasserstein, P)
    e=distanceSqr(Fisher, P)
    metric=Metric(Int(logdet0)) # or metric=logdet0
    s=string(metric) # check what is the current metric
    f=distance²(metric, P) #using the alias distance²

 ## Examples (2)
    using PosDefManifold
    P=randP(10)
    Q=randP(10)
    d=distanceSqr(logEuclidean, P, Q)
    e=distance²(Jeffrey, P, Q)

"""
function distanceSqr(metric::Metric, P::ℍ)

    if      metric==Euclidean
    return  sumOfSqr(P-I)

    elseif  metric==invEuclidean
    return  sumOfSqr(inv(P)-I)

    elseif  metric in (logEuclidean, Fisher)
    return  𝚺(log.(eigvals(P)).^2)

    elseif  metric==logdet0
    return  real(logdet((P+I)/2) - logdet(P)/2)

    elseif  metric==ChoEuclidean
    return  sumOfSqr(choL(P)-I)

    elseif  metric==logCholesky
            LP=choL(P)
    return  sumOfSqrTril(LP, -1) + sumOfSqrDiag(𝑓𝑫(log, LP))

    elseif  metric==Jeffrey
    return  tr(P)/2 + tr(inv(P))/2 - size(P, 1)

    elseif  metric==VonNeumann # see squared distance
            𝓵P=log(P)
    return  real(tr(P*𝓵P) - tr(𝓵P)/2)

    elseif  metric==Wasserstein
    return  tr(P) + size(P, 1) - 2*tr(sqrt(P))

    else    @warn("in RiemannianGeometryP.distanceSqr function
             (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end #function

function distanceSqr(metric::Metric, P::ℍ, Q::ℍ)
    if      metric==Euclidean
    return  sumOfSqr(P - Q)

    elseif  metric==invEuclidean
    return  sumOfSqr(inv(P) - inv(Q))

    elseif  metric==logEuclidean
    return  sumOfSqr(log(P) - log(Q))

    elseif  metric==Fisher
    return  𝚺(log.(eigvals(P, Q)).^2)

    elseif  metric==logdet0
    return  real(logdet((P + Q) / 2) - logdet(P * Q)/2)

    elseif  metric==ChoEuclidean
    return  sumOfSqr(choL(P)-choL(Q))

    elseif  metric==logCholesky
            LP = choL(P)
            LQ = choL(Q)
    return  sumOfSqrTril(tril(LP,-1)-tril(LQ,-1), -1) + sumOfSqrDiag(𝑓𝑫(log, LP)-𝑓𝑫(log, LQ))

    elseif  metric==Jeffrey
            n=size(P, 1)  #using formula tr(Q⁻¹P)/2 + tr(P⁻¹Q)/2 -n
    return  real( (𝚺(colProd(inv(Q), P, i, i) for i=1:n)/2 + 𝚺(colProd(inv(P), Q, i, i) for i=1:n))/2  ) - n
    #return  real(tr(inv(Q)*P)/2 + tr(inv(P)*Q)/2) - size(P, 1)

    elseif  metric==VonNeumann      # using formula: tr(PlogP - PlogQ + QlogQ - QlogP)/2=
            n=size(P, 1)            # (tr(P(logP - LoqQ)) + tr(Q(logQ - logP)))/2=
            R=log(P)-log(Q)         # (tr(P(logP - LoqQ)) - tr(Q(logP - LoqQ)))/2
    return  real( 𝚺(colProd(P, R, i, i) for i=1:n) - 𝚺(colProd(Q, R, i, i) for i=1:n) /2)
    #return  (tr(P*R) - tr(Q*R))/2

    elseif  metric==Wasserstein
            P½=sqrt(P)
    return  tr(P) + tr(Q) -2*real(tr(sqrt(ℍ(P½*Q*P½))))

    else    @warn("in RiemannianGeometryP.distanceSqr function
                    (PosDefManifold Package): the chosen 'metric' does not exist")
    end #if
end # function
distance²=distanceSqr # alias


"""
    (1) distance(metric::Metric, P::ℍ)
    (2) distance(metric::Metric, P::ℍ, Q::ℍ)

 (1) Return ``δ(P, I)``, the *distance* between positive definite matrix ``P`` and
 the identity matrix.

 (2) Return ``δ(P, Q)``, the *distance* between positive definite
 matrices ``P`` and ``Q``.

 This is the square root of [`distanceSqr`](@ref)
 and is invoked with the same syntax therein.

 **See also**: [`distanceMat`](@ref).
"""
distance(metric::Metric, P::ℍ) = √(distanceSqr(metric, P))
distance(metric::Metric, P::ℍ, Q::ℍ) = √(distanceSqr(metric, P, Q))



# -----------------------------------------------------------
# 3. Inter-distance matrix, Laplacian and Spectral Embedding
# -----------------------------------------------------------

# Internal Function for fast computation of inter_distance matrices
function GetdistSqrMat(metric::Metric, 𝐏::ℍVector)
    n, k=_attributes(𝐏)
    △=zeros(k,  k)

    if      metric==invEuclidean
            𝐏𝓲=[inv(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(𝐏𝓲[i] - 𝐏𝓲[j])  end

    elseif  metric==logEuclidean
            𝐏𝓵=[log(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(𝐏𝓵[i] - 𝐏𝓵[j])  end

    elseif  metric==ChoEuclidean
            𝐏L=[choL(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(𝐏L[i] - 𝐏L[j])  end

    elseif  metric==logCholesky
            𝐏L=[choL(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqrTril(tril(𝐏L[i], -1)-tril(𝐏L[j], -1), -1)
                        + sumOfSqrDiag(𝑓𝑫(log, 𝐏L[i])-𝑓𝑫(log, 𝐏L[j])) end

    elseif  metric==Jeffrey
            𝐏𝓲=[inv(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k # optimize computingonly diagonal elements
                #△[i, j]=0.5*(tr(𝐏𝓲[j]*𝐏[i]) + tr(𝐏𝓲[i]*𝐏[j])) - n   end
                △[i, j]=real((𝚺(colProd(𝐏𝓲[j], 𝐏[i], l, l) for l=1:n) + 𝚺(colProd(𝐏𝓲[i], 𝐏[j], l, l) for l=1:n))/2) - n end

    elseif  metric==VonNeumann  # using formula: tr( PlogP + QLoqQ - PlogQ - QlogP)
            𝓵𝐏=[log(P)  for P in 𝐏]
            ℒ=[P*log(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                #△[i, j]=(tr(ℒ[i])+tr(ℒ[j])-tr(𝐏[i] * 𝓵𝐏[j])-tr(𝐏[j] * 𝓵𝐏[i]))/2   end
                △[i, j]=(tr(ℒ[i])+tr(ℒ[j])-real(𝚺(colProd(𝐏[i], 𝓵𝐏[j], l, l) for l=1:n)+𝚺(colProd(𝐏[j], 𝓵𝐏[i], l, l) for l=1:n)))/2 end

    elseif  metric==Wasserstein
            𝐏½=[sqrt(P) for P in 𝐏]
            for j in 1:k-1, i in j+1:k
                △[i, j]=tr(𝐏[i]) + tr(𝐏[j]) -2*real(tr(sqrt(ℍ(𝐏½[i] * 𝐏[j] * 𝐏½[i])))) end

    elseif  metric in (Euclidean, Fisher, logdet0)
            for j in 1:k-1, i in j+1:k
                △[i, j]=distanceSqr(metric, 𝐏[i], 𝐏[j])  end

    else    @warn("in RiemannianGeometryP.distanceSqrMat or .distanceMat function
                         (PosDefManifold Package): the chosen 'metric' does not exist")

    end # If

    return △
end #function



"""
    distanceSqrMat(metric::Metric, 𝐏::ℍVector)

 **alias**: `distance²Mat`

 Given a 1d array `𝐏` of ``k`` positive definite matrices
 ``{P_1,...,P_k}`` of [ℍVector type](@ref), create the ``k⋅k`` real `Hermitian`
 matrix comprising elements ``δ^2(P_i, P_j)\\textrm{, for all }i≠j``.

 This is the matrix of all *squared inter-distances* (zero on diagonal), using the
 specified `metric`, of type [Metric::Enumerated type](@ref),
 giving rise to distance function ``δ``. See [`distanceSqr`](@ref).

 **See**: [distance](@ref).

 **See also**: [`laplacian`](@ref), [`laplacianEigenMaps`](@ref), [`spectralEmbedding`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
    # Compute the squared inter-distance matrix according to the log Euclidean metric.
    # This is much faster as compared to the Fisher metric and in general
    # it is a good approximation.
    Dsqr=distanceSqrMat(logEuclidean, Pset)
    # or, using unicode: Δ²=distanceSqrMat(logEuclidean, 𝐏)

"""
distanceSqrMat(metric::Metric, 𝐏::ℍVector)=ℍ(GetdistSqrMat(metric, 𝐏), :L)
distance²Mat=distanceSqrMat


"""
    distanceMat(metric::Metric, 𝐏::ℍVector)

 Given a 1d array `𝐏` of ``k`` positive definite matrices
 ``{P_1,...,P_k}`` of [ℍVector type](@ref), create the ``k⋅k`` real `Hermitian`
 matrix comprising elements
 ``δ(P_i, P_j)\\textrm{, for all }i≠j``.

 This is the matrix of all *inter-distances* (zero on diagonal), using the
 specified `metric`, of type [Metric::Enumerated type](@ref),
 giving rise to distance ``δ``. See [`distance`](@ref).

 The elements of this matrix are the square root of
 [`distanceSqrMat`](@ref).

 **See**: [distance](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
    D=distanceMat(Fisher, Pset)
    # or, using unicode: Δ=distanceMat(Fisher, 𝐏)
"""
distanceMat(metric::Metric, 𝐏::ℍVector)=ℍ(sqrt.(GetdistSqrMat(metric, 𝐏)), :L)


"""
    laplacian(Δ²)

 Given a matrix of squared inter-distances ``Δ^2``,
 computed for examples by function [`distanceSqrMat`](@ref),
 return the *normalized Laplacian*.

 First, a [Gaussian radial basis functions](https://bit.ly/1HVyf55)
 is applied to all elements of ``Δ^2``, such as

 ``W_{ij} = exp(\\frac{\\displaystyle{-Δ^2_{ij}}}{\\displaystyle{ε}})``,

  where ``ε`` is the Gaussian scale parameter chosen automatically
  as the median of the elements ``Δ^2_{ij}``.

  Finally, the normalized Laplacian is defined as

 ``Ω = D^{-1/2}WD^{-1/2}``,

  where ``D`` is the diagonal matrix holding on the main diagonal
  the sum of the rows (or columns) of ``W``.

!!! note "Nota Bene"
    The normalized Laplacian as here defined can be requested for any
    input matrix of squared inter-distances, for example,
    those obtained on scalars or on vectors using appropriate metrics.

 **See also**: [`distanceSqrMat`](@ref), [`laplacianEigenMaps`](@ref), [`spectralEmbedding`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
    Dsqr=distanceSqrMat(Fisher, Pset) # or: Δ²=distanceSqrMat(Fisher, 𝐏)
    lap=laplacian(Dsqr) # or: Ω=laplacian(Δ²)

 """
function laplacian(Δ²)
    (r, c)=size(Δ²)
    epsilon=median([Δ²[i, j] for j=1:c-1 for i=j+1:r]) # use geometric mean instead
    L=Matrix{eltype(Δ²)}(undef, r, c)
    for i=1:r L[i, i]=1.0 end
    for j=1:c-1, i=j+1:r L[i, j]=exp(-Δ²[i, j]/epsilon)  end
    W=ℍ(L, :L)
    Dnorms=⋱([1/(√(𝚺(W[:, j]))) for j=1:c])
    return ℍ(Dnorms * W * Dnorms) # Ω, see laplacianEigenMaps
end


"""
    laplacianEigenMaps(Ω, q::Int;
                      <tol=1e-9, maxiter=300, ⍰=false>)

 **alias**: `laplacianEM`

 Given a normalized Laplacian ``Ω`` (see [`laplacian`](@ref) ) return
 the *eigen maps* in ``q`` dimensions, i.e., the ``q`` eigenvectors of
 the normalized Laplacian associated with the largest ``q``
 eigenvalues, excluding the first (which is always equal to 1.0).

 The eigenvectors of the normalized Laplacian are computed by the
 power iterations+modified Gram-Schmidt method,
 allowing calling this function even for big Laplacian matrices.

 Return the 4-tuple ``(Λ, U, iterations, convergence)``, where:
 - ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps,
 - ``U`` holds in columns the eigen maps, that is, the ``q`` eigenvectors,
 - ``iterations`` is the number of iterations executed by the power method,
 - ``convergence`` is the convergence attained by the power method.

 The eigenvectors of ``U`` holds the coordinates of the points in a
 low-dimension Euclidean space (typically two or three).
 This is done for, among other purposes, classifying them and
 following their trajectories over time or other dimensions.
 For examples of applications see Ridrigues et *al.* (2018) [🎓](@ref)
 and references therein.

!!! note "Nota Bene"
    The maximum value of ``q`` that can be requested is ``n-1``,
    where ``n`` is the size of the Laplacian.
    In general, ``q=2`` or ``q=3`` is requested.


 **Arguments**: `(Ω, q; <tol=1e-9, maxiter=300, ⍰=false>)`:
 - ``Ω`` is a normalized Laplacian obtained by the [`laplacian`](@ref) function,
 - ``q`` is the dimension of the Laplacian eigen maps;
 - The following are *<optional keyword arguments>* for the power method iterative algorithm:
   * `tol` is the tolerance for convergence,
   * `maxiter` is the maximum number of iterations allowed,
   * if `⍰` is true, the convergence at all iterations will be printed.

 **See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`spectralEmbedding`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
    Dsqr=distanceSqrMat(Fisher, Pset) #or: Δ²=distanceSqrMat(Fisher, 𝐏)
    lap= laplacian(Dsqr) # or: Ω=laplacian(Δ²)
    evalues, maps, iterations, convergence=laplacianEM(lap, 2)
    evalues, maps, iterations, convergence=laplacianEM(lap, 2; maxiter=500)
    evalues, maps, iterations, convergence=laplacianEM(lap, 2; ⍰=true)

"""
function laplacianEigenMaps(Ω, q::Int;
                            tol=1e-9, maxiter=300, ⍰=false)
    (Λ, U, iter, conv) =
        powIter(Ω, q+1; evalues=true, tol=tol, maxiter=maxiter, ⍰=⍰)
    return ⋱(Λ[2:q+1, 2:q+1]), U[1:size(U, 1), 2:q+1], iter, conv
end
laplacianEM=laplacianEigenMaps


"""
    spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int;
                     <tol=1e-9, maxiter=300, ⍰=false>)

 Given a 1d array `𝐏` of ``k`` positive definite matrices ``{P_1,...,P_k}``,
 compute its *eigen maps* in ``q`` dimensions.

 This function runs one after the other the functions:
 - [`distanceSqrMat`](@ref) (compute the squared inter-distance matrix),
 - [`laplacian`](@ref) (compute the normalized Laplacian),
 - [`laplacianEigenMaps`](@ref) (get the eigen maps).

  Return the 4-tuple `(Λ, U, iterations, convergence)`, where:
 - ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps,
 - ``U`` holds in columns the ``q`` eigenvectors, i.e., the ``q`` coordinates of the points in the embedded space,
 - ``iterations`` is the number of iterations executed by the power method,
 - ``convergence`` is the convergence attained by the power method.

 **Arguments** `(metric, 𝐏, q, <tol=1e-9, maxiter=300, ⍰=false>)`:
 - `metric` is the metric of type [Metric::Enumerated type](@ref) used for computing the inter-distances,
 - `𝐏` is a 1d array of ``k`` positive matrices of [ℍVector type](@ref),
 - ``q`` is the dimension of the Laplacian eigen maps;
 - The following are *<optional keyword arguments>* for the power method iterative algorithm:
   * `tol` is the tolerance for convergence of the power method,
   * `maxiter` is the maximum number of iterations allowed for the power method,
   * if `⍰` is true the convergence at all iterations will be printed.

 **See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`laplacianEigenMaps`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    Pset=randP(10, 4) # or, using unicode: 𝐏=randP(10, 4)
    evalues, maps, iterations, convergence=spectralEmbedding(logEuclidean, Pset, 2)
    evalues, maps, iterations, convergence=spectralEmbedding(logEuclidean, Pset, 2; ⍰=true)

"""
function spectralEmbedding(metric::Metric, 𝐏::ℍVector, q::Int;
                           tol=1e-9, maxiter=300, ⍰=false)
    return (Λ, U, iter, conv) =
      laplacianEM(laplacian(distance²Mat(metric, 𝐏)), q; tol=tol, maxiter=maxiter, ⍰=⍰)
end



# -----------------------------------------------------------
# 4. Means (centers of mass, barycenters, ...)
# -----------------------------------------------------------

"""
    (1) mean(metric::Metric, P::ℍ, Q::ℍ)

    (2) mean(metric::Metric, 𝐏::ℍVector;
            <w::Vector=[], ✓w=true>)

 (1) Mean of two positive definite matrices, passed in arbitrary order as
 arguments ``P`` and ``Q``, using the specified `metric` of type
 [Metric::Enumerated type](@ref).
 The order is arbitrary as all metrics implemented in **PosDefManifold** are symmetric.
 This is the midpoint of the geodesic.
 For the weighted mean of two positive definite matrices use instead
 the [`geodesic`](@ref) function.
 ``P`` and ``Q`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 (2) [Fréchet mean](@ref) of an 1d array ``𝐏`` of ``k`` positive definite
 matrices ``𝐏={P_1,...,P_k}`` of [ℍVector type](@ref),
 with optional non-negative real weights ``w={w_1,...,w_k}`` and using the
 specified `metric`as in (1).

 If you don't pass a weight vector with *<optional keyword argument>* ``w``,
 return the *unweighted mean*.

 If *<optional keword argument>* `✓w=true` (default), the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed
 and should be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

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

 and for those that verify an equation:

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

 ## Examples
    using LinearAlgebra, Statistics, PosDefManifold
    # Generate 2 random 3x3 SPD matrices
    P=randP(3)
    Q=randP(3)
    M=mean(logdet0, P, Q) # (1)
    M=mean(logdet0, P, Q) # (1)

    R=randP(3)
    # passing several matrices and associated weights listing them
    # weights vector, does not need to be normalized
    mean(Fisher, ℍVector([P, Q, R]); w=[1, 2, 3])

    # Generate a set of 4 random 3x3 SPD matrices
    Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)
    weights=[1, 2, 3, 1]
    # passing a vector of Hermitian matrices (ℍVector type)
    M=mean(Euclidean, Pset; w=weights) # (2) weighted Euclidean mean
    M=mean(Wasserstein, Pset)  # (2) unweighted Wassertein mean
    # using unicode: M=mean(Wasserstein, 𝐏)
"""
mean(metric::Metric, P::ℍ, Q::ℍ) = geodesic(metric, P, Q, 0.5)

function mean(metric::Metric, 𝐏::ℍVector;
              w::Vector=[], ✓w=true)
    # iterative solutions
    if      metric == Fisher
            (G, iter, conv)=powerMean(𝐏, 0; w=w, ✓w=✓w)
            return G
    elseif  metric == logdet0
            (G, iter, conv)=logdet0Mean(𝐏; w=w, ✓w=✓w)
            return G
    elseif  metric == Wasserstein
            (G, iter, conv)=wasMean(𝐏; w=w, ✓w=✓w)
            return G
    end

    # closed-form expressions
    n, k = _attributes(𝐏)
    isempty(w) ? nothing : v = _getWeights(w, ✓w, k)

    if  metric == Euclidean
        if isempty(w)   return ℍ(𝛍(𝐏))
        else            return ℍ(𝚺(ω*P for (ω, P) in zip(v, 𝐏)))
        end

    elseif metric == invEuclidean
        if isempty(w)   return inv(ℍ(𝛍(inv, 𝐏)))
        else            return inv(ℍ(𝚺(ω*inv(P) for (ω, P) in zip(v, 𝐏))))
        end

    elseif metric == logEuclidean
        if isempty(w)   return ℍ(exp(ℍ(𝛍(log, 𝐏))))
        else            return ℍ(exp(ℍ(𝚺(ω*log(P) for (ω, P) in zip(v, 𝐏)))))
        end

    elseif metric == ChoEuclidean
        if isempty(w)   L = 𝛍(choL, 𝐏)
        else            L = 𝚺(ω*choL(P) for (ω, P) in zip(v, 𝐏))
        end
        return ℍ(L*L')

    elseif metric == logCholesky # Aggiusta!
        L𝐏=[choL(P) for P in 𝐏]
        if isempty(w)
            T=𝛍(tril(L,-1) for L in L𝐏) + exp(mean(𝑓𝑫(log, L) for L in L𝐏))
        else
            T=𝚺(ω*tril(L,-1) for (ω, L) in zip(v, L𝐏))
                + exp(𝚺(ω*𝑓𝑫(log, L) for (ω, L) in zip(v, L𝐏)))
        end
        return ℍ(T*T')

    elseif metric == Jeffrey
        P=mean(Euclidean, 𝐏; w=w, ✓w=✓w)
        Q=mean(invEuclidean, 𝐏; w=w, ✓w=✓w)
        P½, P⁻½=pow(P, 0.5, -0.5)
        return ℍ(P½ * sqrt(ℍ(P⁻½ * Q * P⁻½)) * P½)

    elseif metric == VonNeumann
        @warn "function RiemannianGeometryP.mean and .geodesic not defined for metric $metric"

    else
        @warn "in RiemannianGeometryP.mean function: the chosen 'metric' does not exist"
    end # if metric
end # function

"""
    means(metric::Metric, ℘::ℍVector₂)

 Given a 2d array `℘` of positive definite matrices as an [ℍVector₂ type](@ref)
 compute the [Fréchet mean](@ref) for as many [ℍVector type](@ref) object
 as hold in `℘`, using the specified `metric` of type
 [Metric::Enumerated type](@ref).
  Return the means in a vector of Hermitian matrices, that is, as an `ℍVector` type.

 The weigted Fréchet mean is not supported in this function.

  **See also**: [`mean`](@ref).

  ## Examples
     using PosDefManifold
     # Generate a set of 4 random 3x3 SPD matrices
     Pset=randP(3, 4) # or, using unicode: 𝐏=randP(3, 4)
     # Generate a set of 40 random 4x4 SPD matrices
     Qset=randP(3, 4) # or, using unicode: 𝐐=randP(3, 4)
     # listing directly ℍVector objects
     means(logEuclidean, ℍVector₂([Pset, Qset])) # or: means(logEuclidean, ℍVector₂([𝐏, 𝐐]))
     # note that [𝐏, 𝐐] is actually a ℍVector₂ type object

     # creating and passing an object of ℍVector₂ type
     sets=ℍVector₂(undef, 2) # or: ℘=ℍVector₂(undef, 2)
     sets[1]=Pset # or: ℘[1]=𝐏
     sets[2]=Qset # or: ℘[2]=𝐐
     means(logEuclidean, sets) # or: means(logEuclidean, ℘)

"""
means(metric::Metric, ℘::ℍVector₂)=ℍVector([mean(metric, 𝐏) for 𝐏 in ℘])


"""
    generalizedMean(𝐏::ℍVector, p::Real;
                   <w::Vector=[], ✓w=true>)

 Given a 1d array `𝐏` of ``k`` positive definite matrices ``𝐏={P_1,...,P_k}``
 of [ℍVector type](@ref) and optional non-negative real weights vector ``w={w_1,...,w_k}``,
 return the *weighted generalized mean* ``G`` with real parameter ``p``, that is,

 ``G=\\big(\\sum_{i=1}^{k}w_iP_i^p\\big)^{1/p}``.

 If you don't pass a weight vector with *<optional keyword argument>* ``w``,
 return the *unweighted generalized mean*.

 ``G=\\big(\\sum_{i=1}^{k}P_i^p\\big)^{1/p}``.

 If *<optional keword argument>* `✓w=true` (default), the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed.
 This option is provided to allow
 calling this function repeatedly without normalizing the weights each time.

 The following special cases for parameter ``p`` are noteworthy:
 - For ``p=\\frac{1}{2}`` the generalized mean is the [modified Bhattacharyya mean](@ref).
 - For ``p=1`` the generalized mean is the [Euclidean](@ref) mean.
 - For ``p=-1`` the generalized mean is the [inverse Euclidean](@ref) mean.
 - For ``p=0`` the generalized mean is the [log Euclidean](@ref) mean, which is the [Fisher](@ref) mean when matrices in 𝐏 all pair-wise commute.

 Notice that when matrices in 𝐏 all pair-wise commute,
 the generalized means coincide with the [power means](@ref)
 for any ``p∈[-1, 1]`` and for ``p=0.5`` it coincides also with the
 *Wasserstein* mean. For this reason the generalized means are used
 as default initialization of both the [`powerMean`](@ref) and [`wasMean`](@ref)
 algorithm.

 **See**: [generalized means](@ref).

 **See also**: [`powerMean`](@ref).

 ## Examples
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

"""
function generalizedMean(𝐏::ℍVector, p::Real;
                         w::Vector=[], ✓w=true)
    if     p == -1 return mean(invEuclidean, 𝐏; w=w, ✓w=✓w)
    elseif p ==  0 return mean(logEuclidean, 𝐏; w=w, ✓w=✓w)
    elseif p ==  1 return mean(Euclidean, 𝐏;    w=w, ✓w=✓w)
    else
        n, k=_attributes(𝐏)
        if isempty(w)
            return ℍ(𝛍(P^p for P in 𝐏))^(1/p)
        else
            v=_getWeights(w, ✓w, k)
            return ℍ(𝚺(ω*P^p for (ω, P) in zip(v, 𝐏)))^(1/p)
        end # if w
    end # if p
end # function


"""

    logdet0Mean(𝐏::ℍVector;
               <w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false>)

 Given a 1d array ``𝐏`` of ``k`` positive definite matrices ``𝐏={P_1,...,P_k}``
 of [ℍVector type](@ref) and optional non-negative real weights vector ``w={w_1,...,w_k}``,
 return the 3-tuple ``(G, iter, conv)``, where ``G`` is the mean according
 to the [logdet zero](@ref) metric and ``iter``, ``conv`` are the number of iterations
 and convergence attained by the algorithm.
 Mean ``G`` is the unique positive definite matrix satisfying

 ``\\sum_{i=1}^{k}w_i\\big(\\frac{1}{2}P_i+\\frac{1}{2}G\\big)^{-1}=G^{-1}``.

 For estimating it, this function implements the fixed-point iteration algorithm
suggested by (Moakher, 2012, p315)[🎓](@ref), yielding iterations

 ``G ← \\frac{1}{2}\\big(\\sum_{i=1}^{k}w_i(P_i+G)^{-1}\\big)^{-1}``.

 If you don't pass a weight vector with *<optional keyword argument>* ``w``,
 return the *unweighted logdet zero mean*.

 If *<optional keword argument>* `✓w=true` (default), the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed
 and should be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

 The following are more *<optional keyword arguments*>:
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the [log Euclidean](@ref) mean will be used,
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the saddle point,
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges a **warning** is printed indicating the iteration
    when this happened.

 **See**: [logdet zero](@ref) metric, [modified Bhattacharyya mean](@ref).

 ## Examples
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
    G, iter, conv = logdet0Mean(Pset; w=weights, ⍰=true)

    # now suppose Pset has changed a bit, initialize with G to hasten convergence
    Pset[1]=ℍ(Pset[1]+(randP(3)/100))
    G, iter, conv = logdet0Mean(Pset; w=weights, ✓w=false, ⍰=true, init=G)

"""
function logdet0Mean(𝐏::ℍVector;
                     w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false)
    maxIter=500
    n, k = _attributes(𝐏)
    l=k/2
    isempty(w) ? v=[] : v = _getWeights(w, ✓w, k)
    init == nothing ? M = mean(logEuclidean, 𝐏; w=w, ✓w=false) : M = ℍ(init)
    💡 = similar(M, eltype(M))
    iter, conv, oldconv = 1, 0., maxpos
    ⍰ && @info("Iterating RlogDetMean Fixed-Point...")

    @inbounds while true
        if isempty(w)
            💡 = l * inv(ℍ(𝚺(inv(ℍ(P+M)) for P in 𝐏)))
        else
            💡 = 0.5 * inv(ℍ(𝚺(ω * inv(ℍ(P+M)) for (ω, P) in zip(v, 𝐏))))
        end
        conv = norm(💡-M)/norm(M)
        ⍰ && println("iteration: ", iter, "; convergence: ", conv)
        diverging = conv > oldconv
        diverging ? @warn("logdet0Mean diverged at:", iter) : oldconv=conv
        iter==maxIter || conv <= tol ? break : M = 💡 # diverging ||
        iter += 1
    end # while

    return (💡, iter, conv)
end


"""
    wasMean(𝐏::ℍVector;
           <w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false>)

 Given a 1d array `𝐏` of ``k`` positive definite matrices ``𝐏={P_1,...,P_k}``
 of [ℍVector type](@ref) and optional non-negative real weights vector ``w={w_1,...,w_k}``,
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

 If *<optional keword argument>* `✓w=true` (default), the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed
 and they should be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

 The following are more *<optional keyword arguments*>:
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with ``p=0.5`` will be used,
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the true solution,
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges a **warning** is printed indicating the iteration
    when this happened.

 **See**: [Wasserstein](@ref) metric.

 ## Examples
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
    G, iter, conv = wasMean(Pset; w=weights, ⍰=true)

    # now suppose 𝐏 has changed a bit, initialize with G to hasten convergence
    Pset[1]=ℍ(Pset[1]+(randP(3)/100))
    G, iter, conv = wasMean(Pset; w=weights, ⍰=true, init=G)

"""
function wasMean(𝐏::ℍVector;
                 w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false)

    maxIter=500
    iter, conv, oldconv, maxIter, (n, k) = 1, 0., maxpos, 500, _attributes(𝐏)
    isempty(w) ? v=[] : v = _getWeights(w, ✓w, k)
    init == nothing ? M = generalizedMean(𝐏, 0.5; w=v, ✓w=false) : M = ℍ(init)
    💡 = similar(M, eltype(M))
    ⍰ && @info("Iterating wasMean Fixed-Point...")

    @inbounds while true
        S, W=pow(M, 0.5, -0.5)
        if isempty(w)
            💡 = ℍ(W * sqr(ℍ(𝛍(sqrt(ℍ(S*P*S)) for P in 𝐏))) * W)
        else
            💡 = ℍ(W * sqr(ℍ(𝚺((sqrt(ℍ(S*P*S)) * ω) for (ω, P) in zip(v, 𝐏)))) * W)
        end
        conv = norm(💡-M)/norm(M)
        ⍰ &&  println("iteration: ", iter, "; convergence: ", conv)
        diverging = conv > oldconv
        diverging ? @warn("wasMean diverged at:", iter) : oldconv=conv
        iter==maxIter || conv <= tol ? break : M = 💡 # diverging ||
        iter += 1
    end # while

    return (💡, iter, conv)
end


"""
    powerMean(𝐏::ℍVector, p::Real;
             <w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false>)

 Given a 1d array `𝐏` of ``k`` positive definite matrices ``𝐏={P_1,...,P_k}``
 of [ℍVector type](@ref),
 an optional non-negative real weights vector ``w={w_1,...,w_k}`` and
 a real parameter `p` ``\\in[-1, 1]``, return the
 3-tuple ``(G, iter, conv)``, where ``G`` is
 Lim and Palfia (2012)'s [power means](@ref)  of order ``p`` and
 ``iter``, ``conv`` are the number of iterations
 and convergence attained by the algorithm, respectively.
 Mean ``G`` is the unique positive definite matrix satisfying

 ``G=\\sum_{i=1}^{k}(w_iG\\textrm{#}_pP_i)``,

 where ``G\\textrm{#}_pP_i`` is the [Fisher](@ref) geodesic equation.
 In particular:

 - with ``p=-1`` this is the *harmonic mean* (see the [inverse Euclidean](@ref))
 - with ``p=+1`` this is the *arithmetic mean* (see the [Euclidean](@ref))
 - at the limit of ``p`` evaluated at zero from both side this is the *geometric mean* (see the [Fisher](@ref) metric).

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

 If *<optional keword argument>* `✓w=true` (default), the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed
 and should type be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

 The following are more *<optional keyword arguments*>:
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with parameter ``p`` will be used.
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the true solution;
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges a **warning** is printed indicating the iteration
    when this happened.

 **See**: [power means](@ref), [generalized means](@ref), [modified Bhattacharyya mean](@ref).

 ## Examples
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
    G, iter, conv = powerMean(Pset, 0.5; w=weights, ⍰=true)

    # now suppose 𝐏 has changed a bit, initialize with G to hasten convergence
    Pset[1]=ℍ(Pset[1]+(randP(3)/100))
    G, iter, conv = powerMean(Pset, 0.5; w=weights, ⍰=true, init=G)

"""
function powerMean(𝐏::ℍVector, p::Real;
                   w::Vector=[], ✓w=true, init=nothing, tol=1e-9, ⍰=false)
  if !(-1<=p<=1) @error("The parameter p for power means must be in range [-1...1]")
  else
    if     p ≈-1
            return (mean(invEuclidean, 𝐏; w=w, ✓w=✓w), 1, 0)
    elseif p ≈ 0
            LE=mean(logEuclidean, 𝐏, w=w, ✓w=✓w)
            P, iter1, conv1=powerMean(𝐏,  0.01; w=w, ✓w=✓w, init=LE, tol=tol, ⍰=⍰)
            Q, iter2, conv2=powerMean(𝐏, -0.01; w=w, ✓w=✓w, init=P, tol=tol, ⍰=⍰)
            return (geodesic(Fisher, P, Q,  0.5), iter1+iter2, (conv1+conv2)/2)
    elseif p ≈ 1
                return (mean(Euclidean, 𝐏; w=w, ✓w=✓w), 1, 0)
    else
        # Set Parameters
        n, k = _attributes(𝐏)
        sqrtn, absp, maxIter=√n, abs(p), 500
        r=-0.375/absp
        w≠[] ? v = _getWeights(w, ✓w, k) : v=[]
        init == nothing ? M = generalizedMean(𝐏, p; w=v, ✓w=false) : M = ℍ(init)
        p<0 ? X=ℍ(M^(0.5)) : X=ℍ(M^(-0.5))
        💡, H, 𝒫 = similar(X, eltype(X)), similar(X, eltype(X)), similar(𝐏, eltype(𝐏))
        p<0 ? 𝒫=[inv(P) for P in 𝐏] : 𝒫=𝐏
        iter, conv, oldconv = 1, 0., maxpos
        ⍰ && @info("Iterating powerMean Fixed-Point...")

        @inbounds while true
            if isempty(w)
                H=ℍ(𝛍(pow(ℍ(X*P*X), absp) for P in 𝒫))
            else
                H=ℍ(𝚺(ω * pow(ℍ(X*P*X), absp) for (ω, P) in zip(v, 𝒫)))
            end
            💡=(pow(H, r))*X
            conv=norm(H-I)/sqrtn # relative difference to identity
            ⍰ &&  println("iteration: ", iter, "; convergence: ", conv)
            diverging = conv > oldconv
            diverging ? @warn("powerMean diverged at:", iter) : oldconv=conv
            iter==maxIter || conv <= tol ? break : X = 💡 # diverging ||
            iter += 1
        end # while
    end # if

    if p<0  return ( ℍ((💡)'*💡), iter, conv )
    else    return ( inv(ℍ((💡)'*💡)), iter, conv ) end
  end # if !(-1<=p<=1)
end



# -----------------------------------------------------------
# 5. Tangent Space Tools
# -----------------------------------------------------------

"""
    logMap(metric::Metric, P::ℍ, G::ℍ)

 *Logaritmic Map:* map a positive definite matrix ``P`` from the SPD or
 Hermitian manifold into the tangent space at base-point ``G`` using the [Fisher](@ref) metric.

 ``P`` and ``G`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 The map is defined as

 `` Log_G(P)=S=G^{1/2}\\textrm{log}\\big(G^{-1/2}PG^{-1/2}\\big)G^{1/2}``.

 The result is an `Hermitian` matrix.
 The inverse operation is [`expMap`](@ref).

 **Arguments** `(metric, P, G)`:
 - `metric` is a metric of type [Metric::Enumerated type](@ref).
 - ``P`` is the positive definite matrix to be projected onto the tangent space,
 - ``G`` is the tangent space base point,

 Currently only the [Fisher](@ref) metric is supported for tangent space operations.

 **See also**: [`vecP`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(3)
    Q=randP(3)
    metric=Fisher
    G=mean(metric, P, Q)
    # projecting P at the base point given by the geometric mean of P and Q
    S=logMap(metric, P, G)
"""
function logMap(metric::Metric, P::ℍ, G::ℍ)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * log(ℍ(G⁻½ * P * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.logMap function:
                 only the Fisher metric is supported for the logarithmic map."
    end
end

"""

    expMap(metric::Metric, S::ℍ, G::ℍ)

 *Exponential Map:* map an `Hermitian` matrix ``S`` from the tangent space at base
 point ``G`` into the SPD or Hermitian manifold (using the [Fisher](@ref) metric).

 ``S`` and ``G`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 The map is defined as

 `` Exp_G(S)=P=G^{1/2}\\textrm{exp}\\big(G^{-1/2}SG^{-1/2}\\big)G^{1/2}``.

 The result is a positive definite matrix.
 The inverse operation is [`logMap`](@ref).

 **Arguments** `(metric, S, G)`:
 - `metric` is a metric of type [Metric::Enumerated type](@ref),
 - ``S`` is a Hermitian matrix, real or complex, to be projected on the SPD or Hermitian manifold,
 - ``G`` is the tangent space base point.

  Currently only the Fisher metric is supported for tangent space operations.

 ## Examples
    using PosDefManifold, LinearAlgebra
    P=randP(3)
    Q=randP(3)
    G=mean(Fisher, P, Q)
    # projecting P on the tangent space at the Fisher mean base point G
    S=logMap(Fisher, P, G)
    # adding the identity in the tangent space and reprojecting back onto the manifold
    H=expMap(Fisher, ℍ(S+I), G)
"""
function expMap(metric::Metric, S::ℍ, G::ℍ)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * exp(ℍ(G⁻½ * S * G⁻½)) * G½)
    else @warn "in RiemannianGeometryP.expMap function:
              only the Fisher metric is supported for the exponential map"
    end
end


"""
    vecP(S::ℍ)

 *Vectorize* a tangent vector (matrix) ``S`` (*i.e.*, an `Hermitian` matrix):  mat -> vec.

 It gives weight ``1`` to diagonal elements and √2 to off-diagonal elements
 (Barachant et *al.*, 2012)[🎓](@ref).

 The result is a vector holding ``n(n+1)/2`` elements, where ``n``
 is the size of ``S``.

 ``S`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 The inverse operation is provided by [`matP`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(3)
    Q=randP(3)
    G=mean(Fisher, P, Q)
    # projecting P at the base point given by the geometric mean of P and Q
    S=logMap(Fisher, P, G)
    # vectorize S
    v=vecP(S)
"""
vecP(S::ℍ)=[(if i==j return S[i, j] else return (S[i, j])*sqrt2 end) for j=1:size(S, 2) for i=j:size(S, 1)]


"""
    matP(ς::Vector)

 *Matrizize* a tangent vector (vector) ς :  vec -> mat.

 This is the function reversing the [`vecP`](@ref) function,
 thus the weighting applied therein is reversed as well.

 If ``ς=vecP(S)`` and ``S`` is a ``n⋅n`` Hermitian matrix,
 ``ς``  is a tangent vector of size ``n(n+1)/2``.
 The result of calling ``matP(ς)`` is then ``n⋅n`` matrix ``S``.

 P.S.: This function needs to be rewritten more efficiently

 ## Examples
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
"""
function matP(ς::Vector)
  n=Int((-1+√(1+8*length(ς)))/2) # Size of the matrix whose vectorization vector v has size length(v)
  S=Matrix{eltype(ς)}(undef, n, n)
  l=0
  @inbounds for j in 1:n-1
    l=l+1
    S[j, j]=ς[l]
    for i in j+1:n
      l=l+1
      S[i, j]=invsqrt2*ς[l]
    end
  end
  S[n, n]=ς[end]
  return ℍ(S, :L)
end



# -----------------------------------------------------------
# 6. Procrustes Problems
# -----------------------------------------------------------

"""
    procrustes(P::ℍ, Q::ℍ, extremum="min")

 Given two positive definite matrices ``P`` and ``Q``,
 return by default the solution of problem

 ``\\textrm{argmin}_Uδ(P,U^*QU)``,

 where ``U`` varies over the set of unitary matrices ``𝐔`` and ``δ(.,.)`` is a
 distance or divergence function.
 ``U^*QU`` is named in physics the *unitary orbit* of ``Q``.

 If the argument 'extremum' is passed as "max", it returns instead the solution of

 ``\\textrm{argmax}_Uδ(P,U^*QU)``.

  ``P`` and ``Q`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 As it has been shown in Bhatia and Congedo (2019)[🎓](@ref),
 using each of the [Fisher](@ref), [logdet zero](@ref), [Wasserstein](@ref)
 and the Kullback-Leibler divergence (see [logdet α](@ref)),
 the best approximant to ``P`` from the unitary orbit of ``Q``
 commutes with ``P`` and, surprisingly, has the same closed-form expression, namely

 ``U_Q^↓U_P^{↓*}`` for the argmin and ``U_Q^↑U_P^{↓*}`` for the argmax,

 where ``U^↓`` denotes the eigenvector matrix of the subscript argument with
 eigenvectors in columns sorted by *decreasing* order of corresponding eigenvalues and
 ``U^↑`` denotes the eigenvector matrix of the subscript argument with
 eigenvectors in columns sorted by *increasing* order of corresponding eigenvalues.

 The same solutions are known since a long time also by solving the extremal
 problem here above using the [Euclidean](@ref) metric (Umeyama, 1988).

 ## Examples
    using PosDefManifold
    P=randP(3)
    Q=randP(3)
    # argmin problem
    U=procrustes(P, Q)
    # argmax problem
    V=procrustes(P, Q, "max")
"""
function procrustes(P::ℍ, Q::ℍ, extremum="min")
    Pup=eigvecs(P)
    Qup=eigvecs(Q)
    Pdown=reverse(Pup, dims=(2))
    if      extremum=="min"
            Qdown=reverse(Qup, dims=(2))
            return Qdown*Pdown
    elseif  extremum=="max"
            return Qup*Pdown
    else    @warn "in RiemannianGeometryP.procrustes: the argument 'extremum' is incorrect."
    end
end
