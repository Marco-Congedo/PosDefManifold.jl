#    Unit riemannianGeometry.jl, part of PosDefManifold Package for julia language
#    v 0.1.0 - last update 10th of April 2019
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
#    1. Geodesic Equations
#    2. Distances
#    3. Inter-distance matrices, Laplacian and Spectral Embedding
#    4. Means (center of mass, barycenters, ...)
#    5. Tangent Space
#    6. Procrustes Problems

# -----------------------------------------------------------
# 1. Geodesic Equations
# -----------------------------------------------------------

"""
    geodesic(P::ℍ, Q::ℍ, a::Real, metric::Metric=Fisher)

  Move along the [geodesic](@ref) from point ``P`` to point ``Q``
  (two positive definite matrices) with *arclegth* ``0<=a<=1``,
 using the specified metric, of type [Metric::Enumerated type](@ref).
 By default de [Fisher](@ref) metric is adopted.

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
 so in this case the geodesic is found as the weighted mean [`meanP(@ref)`].
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
|Jeffrey | uses weighted mean [`meanP`](@ref) |
|VonNeumann | N.A.|
|Wasserstein| ``b^2P+a^2Q +ab\\big[(PQ)^{1/2} +(QP)^{1/2}\\big]``|

  **legend:** ``L_X``, ``S_X`` and ``D_X``
   are the Cholesky lower triangle of ``X``, its strictly lower triangular part
   and diagonal part, respectively (hence, ``S_X+D_X=L_X``,  ``L_XL_X^*=X``).

 **See also**: [`meanP`](@ref)

 ## Examples
    using PosDefManifold
    P=randP(10);
    Q=randP(10);
    # Wasserstein mean
    M=geodesic(P, Q, 0.5, Wasserstein)
    # extrapolate suing the Fisher metric
    E=geodesic(P, Q, 2)

"""
function geodesic(P::ℍ, Q::ℍ, a::Real, metric::Metric=Fisher)
    if a ≈ 0 return P end
    if a ≈ 1 return Q end
    b = 1-a

    if      metric==Euclidean
    return  ℍ(P*b + Q*a)

    elseif  metric==invEuclidean
    return  ℍ( inv( ℍ(inv(P)*b + inv(Q)*a) ) )

    elseif  metric==logEuclidean
    return  ℍ( exp( ℍ(log(P)*b + log(Q)*a) ) )

    elseif  metric==Fisher
            P½, P⁻½ = pow(P, 0.5, -0.5)
    return  ℍ( P½ * ℍ((P⁻½ * Q * P⁻½)^a) * P½ )

    elseif  metric in (logdet0, Jeffrey)
    return  meanP([P, Q], metric, w=[b, a], ✓w=false) #! 2

    elseif  metric==VonNeumann
            @warn("An expression for the geodesic is not available for the Von neumann metric")

    elseif  metric==ChoEuclidean
            L=choL(P)*b + choL(Q)*a
    return  ℍ(L*L')

    elseif  metric==logCholesky
            LP=choL(P); LQ=choL(Q); slLP=tril(LP,-1)
            #L=slLP+a*(tril(LQ,-1)-slLP)+⋱(LP)*exp(a*(𝑓𝑫(LQ, log)-𝑓𝑫(LP, log)))
            L=slLP+(tril(LQ,-1)-slLP)*a
            for i in 1:size(P, 1) L[i, i]+=LP[i, i]*exp( (log(LQ[i, i])-log(LP[i, i]))*a ) end
    return  ℍ(L*L')

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
    (1) distanceSqr(P::ℍ, metric::Metric=Fisher
    (2) distanceSqr(P::ℍ, Q::ℍ, metric::Metric=Fisher)

 **alias**: `distance²`

 (1) Return ``δ^2(P, I)``, the *square of the distance* (or *divergence*) of positive definite
 matrix ``P`` from the the identity matrix. See [distance from the origin](@ref).

 (2) Return ``δ^2(P, Q)``, the *square of the distance* (or *divergence*) between two
 positive definite matrices ``P`` and ``Q``. See [distance](@ref).

 In both cases the distance function ``δ`` is induced by the argument `metric` of type
 [Metric::Enumerated type](@ref). By default, the [Fisher](@ref) metric is adopted.

 The distance is always real, non-negative and equal to zero if and only if
 (1) ``P=I`` or (2) ``P=Q``.

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

 **See also**: [`distanceSqrMat`](@ref)

 ## Examples (1)
    using PosDefManifold
    P=randP(10);
    d=distanceSqr(P, Wasserstein)
    e=distanceSqr(P) # uses the default metric (Fisher)
    metric=Metric(Int(logdet0)) # or metric=logdet0
    s=string(metric) # check what is the current metric
    f=distance²(P, metric) #using the alias distance²

 ## Examples (2)
    using PosDefManifold
    P=randP(10);
    Q=randP(10);
    d=distanceSqr(P, Q, Wasserstein)
    e=distance²(P, Q, Jeffrey)

"""
function distanceSqr(P::ℍ, metric::Metric=Fisher)
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
            n=size(P, 1)
    return  sumOfSqrTril(tril(LP,-1), -1) + 𝚺(log(LP[i, i])^2 for i in 1:n)

    elseif  metric==Jeffrey
    return  tr(P)/2 + tr(inv(P))/2 - size(P, 1)

    elseif  metric==VonNeumann # see squared distance
            𝓵P=log(P)
    return  (tr(P*𝓵P) - tr(𝓵P))/2

    elseif  metric==Wasserstein
    return  tr(P+I) - 2*tr(sqrt(P))

    else    @warn("in RiemannianGeometryP.distanceSqr function
             (PosDefManifold Package): the chosen 'metric' does not exist")
    end # if
end #function

function distanceSqr(P::ℍ, Q::ℍ, metric::Metric=Fisher)
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
            n=size(P, 1)
    return  sumOfSqrTril(tril(LP,-1)-tril(LQ,-1), -1)+𝚺((log(LP[i, i])-log(LQ[i, i]))^2 for i in 1:n)

    elseif  metric==Jeffrey
    return  real(tr(inv(Q)*P)/2 + tr(inv(P)*Q)/2) - size(P, 1)

    elseif  metric==VonNeumann      # using formula: tr(PlogP - PlogQ + QlogQ - QlogP)/2=
            𝓵Pm𝓵Q=log(P)-log(Q);         # (tr(P(logP - LoqQ)) + tr(Q(logQ - logP)))/2=
    return  (tr(P*𝓵Pm𝓵Q) - tr(Q*𝓵Pm𝓵Q))/2     # (tr(P(logP - LoqQ)) - tr(Q(logP - LoqQ)))/2

    elseif  metric==Wasserstein
            P½=sqrt(P);
    return  tr(P) + tr(Q) -2tr(sqrt(ℍ(P½*Q*P½)))

    else    @warn("in RiemannianGeometryP.distanceSqr function
                    (PosDefManifold Package): the chosen 'metric' does not exist")
    end #if
end # function
distance²=distanceSqr # alias


"""
    (1) distance(P::ℍ, metric::Metric=Fisher
    (2) distance(P::ℍ, Q::ℍ, metric::Metric=Fisher)

 (1) Return ``δ(P, I)``, the *distance* between positive definite matrix ``P`` and
 the identity matrix.

 (2) Return ``δ(P, Q)``, the *distance* between positive definite
 matrices ``P`` and ``Q``.

 This is the square root of [`distanceSqr`](@ref)
 and is invoked with the same syntax therein.

 **See also**: [`distanceMatrix`](@ref)
"""
distance(P::ℍ, metric::Metric=Fisher) = √(distanceSqr(P, metric))
distance(P::ℍ, Q::ℍ, metric::Metric=Fisher) = √(distanceSqr(P, Q, metric))

# -----------------------------------------------------------
# 3. Inter-distance matrix, Laplacian and Spectral Embedding
# -----------------------------------------------------------

# Internal Function for fast computation of inter_distance matrices
function GetdistanceSqrMat(℘, metric::Metric=Fisher)
    k=length(℘)
    n=size(℘[1], 1)
    △=zeros(k,  k)

    if      metric==invEuclidean
            ℘𝓲=[inv(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(℘𝓲[i] - ℘𝓲[j])  end

    elseif  metric==logEuclidean
            ℘𝓵=[log(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(℘𝓵[i] - ℘𝓵[j])  end

    elseif  metric==ChoEuclidean
            ℘L=[choL(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqr(℘L[i] - ℘L[j])  end

    elseif  metric==logCholesky
            ℘L=[choL(P)     for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=sumOfSqrTril(℘L[i]-℘L[j], -1) + 𝚺((log(℘L[i][l, l])-log(℘L[j][l, l]))^2 for l in 1:n) end

    elseif  metric==Jeffrey
            ℘𝓲=[inv(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=tr(℘𝓲[j]*℘[i])/2 + tr(℘𝓲[i]*℘[j])/2 - n   end

    elseif  metric==VonNeumann  # using formula: tr( PlogP + QLoqQ - PlogQ - QlogP)
            𝓵℘=[log(P)      for P in ℘]
            ℘i𝓵℘i=[P*log(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=(tr(℘i𝓵℘i[i])+tr(℘i𝓵℘i[j])-tr(℘[i] * 𝓵℘[j])-tr(℘[j] * 𝓵℘[i]))/2   end

    elseif  metric==Wasserstein
            ℘½=[sqrt(P) for P in ℘]
            for j in 1:k-1, i in j+1:k
                △[i, j]=tr(℘[i]) + tr(℘[j]) -2*tr(sqrt(℘½[i] * ℘[j] * ℘½[i]'))     end

    elseif  metric in (Euclidean, Fisher, logdet0)
            for j in 1:k-1, i in j+1:k
                △[i, j]=distanceSqr(℘[i], ℘[j], metric)  end

    else    @warn("in RiemannianGeometryP.distanceSqrMat or .distanceMatrix function
                         (PosDefManifold Package): the chosen 'metric' does not exist")

    end # If

    return △
end #function



"""
    distanceSqrMat(℘, metric::Metric=Fisher)

 **alias**: `distance²Mat`

 Given a 1d array `℘` of ``k`` positive definite matrices
 ``{P_1,...,P_k}``, create the ``k⋅k`` real `Hermitian` matrix comprising
 elements ``δ^2(P_i, P_j)\\textrm{, for all }i≠j``.

 This is the matrix of all *squared inter-distances* (zero on diagonal), using the
 specified `metric`, of type [Metric::Enumerated type](@ref),
 giving rise to distance function ``δ``. See [`distanceSqr`](@ref).
 By default, the [Fisher](@ref) metric is adopted.

 **See also**: [`laplacian`](@ref), [`laplacianEigenMaps`](@ref), [`spectralEmbedding`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    ℘=randP(10, 4)
    # Compute the squared inter-distance matrix according to the log Euclidean metric.
    # This is much faster as compared to the Fisher metric and in general
    # it is a good approximation.
    Δ²=distanceSqrMat(℘, logEuclidean)

"""
distanceSqrMat(℘, metric::Metric=Fisher)=ℍ(GetdistanceSqrMat(℘, metric), :L)
distance²Mat=distanceSqrMat


"""
    distanceMatrix(℘, metric::Metric=Fisher)

 **alias**: `distanceMat`

 Given a 1d array `℘` of ``k`` positive definite matrices
 ``{P_1,...,P_k}``, create the ``k⋅k`` real `Hermitian` matrix comprising elements
 ``δ(P_i, P_j)\\textrm{, for all }i≠j``.

 This is the matrix of all *inter-distances* (zero on diagonal), using the
 specified `metric`, of type [Metric::Enumerated type](@ref),
 giving rise to distance ``δ``. See [`distance`](@ref).
 By default, the [Fisher](@ref) metric is adopted.

 The elements of this matrix are the square root of
 [`distanceSqrMat`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    ℘=randP(10, 4)
    Δ=distanceMatrix(℘)
"""
distanceMatrix(℘, metric::Metric=Fisher)=ℍ(sqrt.(GetdistanceSqrMat(℘, metric)), :L)
distanceMat=distanceMatrix


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
    ℘=randP(10, 4)
    Δ²=distanceSqrMat(℘)
    Ω=laplacian(Δ²) # or, equivalently, Ω=RΩ(Δ)

 """
function laplacian(Δ²)
    (r, c)=size(Δ²)
    epsilon=median([Δ²[i, j] for j=1:c-1 for i=j+1:r]) # use geometric mean instead
    L=Matrix{eltype(Δ²)}(undef, r, c)
    for i=1:r L[i, i]=1.0 end
    for j=1:c-1, i=j+1:r L[i, j]=exp(-Δ²[i, j]/epsilon)  end
    W=ℍ(L, :L)
    Dnorms=Diagonal([1/(√(𝚺(W[:, j]))) for j=1:c])
    return ℍ(Dnorms * W * Dnorms) # Ω
end


"""
    laplacianEigenMaps(Ω, q::Int; <tol=1e-9, maxiter=300, ⍰=false>)

 **alias**: `laplacianEM`

 Given a normalized Laplacian ``Ω`` (see [`laplacian`](@ref) ) return
 the *eigen maps* in ``q`` dimensions, i.e., the ``q`` eigenvectors of
 the normalized Laplacian associated with the largest ``q``
 eigenvalues, excluding the first (which is always equal to 1.0).

 The eigenvectors of the normalized Laplacian are computed by the
 power iterations+modified Gram-Schmidt method,
 allowing calling this function even for big Laplacian matrices.

 Return the 4-tuple ``(Λ, U, iterations, convergence)``, where:
 - ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps;
 - ``U`` holds in columns the eigen maps, that is, the ``q`` eigenvectors
 - ``iterations`` is the number of iterations executed by the power method;
 - ``convergence`` is the convergence attained by the power method;

 The eigenvectors of ``U`` holds the coordinates of the points in the
 embedded space. For examples of applications see Ridrigues et al. (2018) [🎓](@ref).

!!! note "Nota Bene"
    The maximum value of ``q`` that can be requested is ``n-1``,
    where ``n`` is the size of the Laplacian.
    In general, ``q=2`` or ``q=3`` is requested.


 **Arguments**: `(Ω, q; <tol=1e-9, maxiter=300, ⍰=false>)`:
 - ``Ω`` is a normalized Laplacian obtained by the [`laplacian`](@ref) function;
 - ``q`` is the dimension of the Laplacian eigen maps;
 - The following are *<optional keyword arguments>* for the power method iterative algorithm:
   * `tol` is the tolerance for convergence;
   * `maxiter` is the maximum number of iterations allowed;
   * if `⍰` is true, the convergence at all iterations will be printed.

 **See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`spectralEmbedding`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    ℘=randP(10, 4)
    evalues, maps, iterations, convergence=laplacianEM(Ω, 2)
    evalues, maps, iterations, convergence=laplacianEM(Ω, 2, maxiter=500)
    evalues, maps, iterations, convergence=laplacianEM(Ω, 2, ⍰=true)

"""
function laplacianEigenMaps(Ω, q::Int; tol=1e-9, maxiter=300, ⍰=false)
    (Λ, U, iter, conv) =
        powIter(Ω, q+1; evalues=true, tol=tol, maxiter=maxiter, ⍰=⍰)
    return Diagonal(Λ[2:q+1, 2:q+1]), U[1:size(U, 1), 2:q+1], iter, conv
end;
laplacianEM=laplacianEigenMaps


"""
    spectralEmbedding(℘, q::Int, metric::Metric=Fisher;
                        <tol=1e-9, maxiter=300, ⍰=false>)

 **alias**: `Rse`

 Given a 1d array `℘` of ``k`` positive definite matrices ``{P_1,...,P_k}``,
 compute its *eigen maps* in ``q`` dimensions.

 This function runs one after the other the functions:
 - [`distanceSqrMat`](@ref) (compute the squared inter-distance matrix),
 - [`laplacian`](@ref) (compute the normalized Laplacian),
 - [`laplacianEigenMaps`](@ref) (get the eigen maps).

  Return the 4-tuple `(Λ, U, iterations, convergence)`, where:
 - ``Λ`` is a ``q⋅q`` diagonal matrix holding on diagonal the eigenvalues corresponding to the ``q`` dimensions of the Laplacian eigen maps;
 - ``U`` holds in columns the ``q`` eigenvectors, i.e., the ``q`` coordinates of the points in the embedded space.
 - ``iterations`` is the number of iterations executed by the power method;
 - ``convergence`` is the convergence attained by the power method;

 **Arguments** `(℘, q, metric, <tol=1e-9, maxiter=300, ⍰=false>)`:
 - `℘` is a 1d array of ``k`` positive matrices;
 - ``q`` is the dimension of the Laplacian eigen maps;
 - `metric` is a metric of type [Metric::Enumerated type](@ref),
   used for computing the inter-distances. By default, the [Fisher](@ref) metric is adopted.
 - The following are *<optional keyword arguments>* for the power method iterative algorithm:
   * `tol` is the tolerance for convergence of the power method;
   * `maxiter` is the maximum number of iterations allowed for the power method;
   * if `⍰` is true the convergence at all iterations will be printed.

 **See also**: [`distanceSqrMat`](@ref), [`laplacian`](@ref), [`laplacianEigenMaps`](@ref).

 ## Examples
    using PosDefManifold
    # Generate a set of 4 random 10x10 SPD matrices
    ℘=randP(10, 4)
    evalues, maps, iterations, convergence=spectralEmbedding(℘, 2)
    evalues, maps, iterations, convergence=spectralEmbedding(℘, 2, ⍰=true)

"""
function spectralEmbedding(℘, q::Int, metric::Metric=Fisher;
                            tol=1e-9, maxiter=300, ⍰=false)
    return (Λ, U, iter, conv) =
            laplacianEM(laplacian(distance²Mat(℘, metric)), q; tol=tol, maxiter=maxiter, ⍰=⍰)
end


# -----------------------------------------------------------
# 4. Means (centers of mass, barycenters, ...)
# -----------------------------------------------------------

# Internal functions
Attributes(℘)=( size(℘[1], 1), length(℘))
#Attributes(℘)=(℘[1] isa Hermitian ? Hermitian : Symmetric, size(℘[1], 1), length(℘))
#Attributes(P::HermOrSym)=(P isa Hermitian ? Hermitian : Symmetric, size(P, 1))


function DoNothing
end

function GetWeights(w::Vector, ✓w::Bool, k::Int)
    if ✓w==true
        s=𝚺(w)
        if s ≉  1.0 return w./s else return w end
    else return w
    end
end

"""
    generalizedMean(℘, p::Real; <w::Vector=[], ✓w::Bool=true>)

 Given a 1d array `℘` of ``k`` positive definite matrices``{P_1,...,P_k}``
 and optional non-negative real weights vector ``w={w_1,...,w_k}``,
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
 - For ``p=0`` the generalized mean is the [log Euclidean](@ref) mean, which is the [Fisher](@ref) mean when matrices in `℘` all pair-wise commute.

 Notice that when matrices in `℘` all pair-wise commute,
 the generalized means coincide with the [power means](@ref)
 for any ``p∈[-1, 1]`` and for ``p=0.5`` it coincides also with the
 *Wasserstein* mean (see [`wasMean`](@ref)). For this reason the generalized means are used
 as default initialization of both the [`powerMean`](@ref) and [`wasMean`](@ref)
 algorithm.

 **See also**: [`powerMean`](@ref)

 ## Examples
    using LinearAlgebra, Statistics, PosDefManifold
    # Generate a set of 4 random 3x3 SPD matrices
    ℘=randP(3, 4)

    # weights vector, does not need to be normalized
    weights=[1, 2, 3, 1]

    # unweighted mean
    G = generalizedMean(℘, 0.25)

    # weighted mean
    G = generalizedMean(℘, 0.5; w=weights)

    # with weights previously normalized we can set ✓w=false
    weights=weights./mean(weights)
    G = generalizedMean(℘, 0.5; w=weights, ✓w=false)

"""
function generalizedMean(℘, p::Real; w::Vector=[], ✓w::Bool=true)
    if     p == -1 return meanP(℘, invEuclidean; w=w, ✓w=✓w)
    elseif p ==  0 return meanP(℘, logEuclidean; w=w, ✓w=✓w)
    elseif p ==  1 return meanP(℘, Euclidean;    w=w, ✓w=✓w)
    else
        n, k=Attributes(℘)
        if isempty(w)
            return ℍ(ℍ(𝛍(P^p for P in ℘))^(1/p))
        else
            v=GetWeights(w, ✓w, k)
            return ℍ(ℍ(𝚺(ω*P^p for (ω, P) in zip(v, ℘)))^(1/p))
        end # if w
    end # if p
end # function


"""

    logdet0Mean(℘; <w::Vector=[], ✓w::Bool=true, init=nothing,
                     tol=1e-9, ⍰=false>)

 Given a 1d array ``℘`` of ``k`` positive definite matrices ``{P_1,...,P_k}``
 and optional non-negative real weights vector ``w={w_1,...,w_k}``,
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
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the [log Euclidean](@ref) mean will be used;
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the saddle point;
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges a **warning** is printed indicating the iteration
    when this happened and the algorithm is interrupted.

 **See**: [logdet zero](@ref) metric, [modified Bhattacharyya mean](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    # Generate a set of 4 random 3x3 SPD matrices
    ℘=randP(3, 4)

    # unweighted mean
    G, iter, conv = logdet0Mean(℘)

    # weights vector, does not need to be normalized
    weights=[1, 2, 3, 1]

    # weighted mean
    G, iter, conv = logdet0Mean(℘, w=weights)

    # print the convergence at all iterations
    G, iter, conv = logdet0Mean(℘, w=weights, ⍰=true)

    # now suppose ℘ has changed a bit, initialize with G to hasten convergence
    ℘[1]=ℍ(℘[1]+(randP(3)/100))
    G, iter, conv = logdet0Mean(℘, w=weights, ✓w=false, ⍰=true, init=G)

"""
function logdet0Mean(℘;    w::Vector=[], ✓w::Bool=true, init=nothing,
                            tol=1e-9, ⍰=false)
    maxIter=500
    n, k = Attributes(℘)
    l=k/2
    isempty(w) ? v=[] : v = GetWeights(w, ✓w, k)
    init == nothing ? M = meanP(℘, logEuclidean, w=w, ✓w=false) : M = ℍ(init)
    M◇ = similar(M, eltype(M))
    iter = 1
    conv = 0.; oldconv=maxpos
    ⍰ && @info("Iterating RlogDetMean Fixed-Point...")

    @inbounds while true
        if isempty(w)
            M◇ = ℍ(l * inv(ℍ(𝚺(inv(ℍ(P+M)) for P in ℘))))
        else
            M◇ = ℍ(0.5 * inv(ℍ(𝚺(ω * inv(ℍ(P+M)) for (ω, P) in zip(v, ℘)))))
        end
        conv = norm(M◇-M)/norm(M)
        ⍰ && println("iteration: ", iter, "; convergence: ", conv)
        diverging = conv > oldconv
        diverging ? @warn("logdet0Mean diverged at:", iter) : oldconv=conv
        iter==maxIter || diverging || conv <= tol ? break : M = M◇
        iter += 1
    end # while

    return (M◇, iter, conv)
end


"""
    wasMean(℘; <w::Vector=[], ✓w::Bool=true, init=nothing,
                 tol=1e-9, ⍰=false>)

 Given a 1d array `℘` of ``k`` positive definite matrices ``{P_1,...,P_k}``
 and optional non-negative real weights vector ``w={w_1,...,w_k}``,
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
 and Metric::Enumerated type be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

 The following are more *<optional keyword arguments*>:
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with ``p=0.5`` will be used;
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the true solution;
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
    In normal circumstances this algorithm converges monothonically.
    If the algorithm diverges a **warning** is printed indicating the iteration
    when this happened and the algorithm is interrupted.

 **See**: [Wasserstein](@ref) metric

 ## Examples
    using LinearAlgebra, PosDefManifold
    # Generate a set of 4 random 3x3 SPD matrices
    ℘=randP(3, 4)

    # unweighted mean
    G, iter, conv = wasMean(℘)

    # weights vector, does not need to be normalized
    weights=[1, 2, 3, 1]

    # weighted mean
    G, iter, conv = wasMean(℘, w=weights)

    # print the convergence at all iterations
    G, iter, conv = wasMean(℘, w=weights, ⍰=true)

    # now suppose ℘ has changed a bit, initialize with G to hasten convergence
    ℘[1]=ℍ(℘[1]+(randP(3)/100))
    G, iter, conv = wasMean(℘, w=weights, ⍰=true, init=G)

"""
function wasMean(℘;    w::Vector=[], ✓w::Bool=true, init=nothing,
                        tol=1e-9, ⍰=false)
    maxIter=500
    n, k = Attributes(℘)
    isempty(w) ? v=[] : v = GetWeights(w, ✓w, k)
    init == nothing ? M = generalizedMean(℘, 0.5; w=v, ✓w=false) : M = ℍ(init)
    M◇ = similar(M, eltype(M))
    iter = 1
    conv = 0.; oldconv=maxpos
    ⍰ && @info("Iterating wasMean Fixed-Point...")

    @inbounds while true
        S, W=pow(M, 0.5, -0.5)
        if isempty(w)
            M◇ = ℍ(W * sqr(ℍ(𝛍(sqrt(ℍ(S*P*S)) for P in ℘))) * W)
        else
            M◇ = ℍ(W * sqr(ℍ(𝚺((sqrt(ℍ(S*P*S)) * ω) for (ω, P) in zip(v, ℘)))) * W)
        end
        conv = norm(M◇-M)/norm(M)
        ⍰ &&  println("iteration: ", iter, "; convergence: ", conv)
        diverging = conv > oldconv
        diverging ? @warn("wasMean diverged at:", iter) : oldconv=conv
        iter==maxIter || diverging || conv <= tol ? break : M = M◇
        iter += 1
    end # while

    return (M◇, iter, conv)
end


"""
    powerMean(℘, p::Real; <w::Vector=[], ✓w::Bool=true, init=nothing,
                            tol=1e-9, ⍰=false>)

 Given a 1d array `℘` of ``k`` positive definite matrices ``{P_1,...,P_k}``,
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
 and Metric::Enumerated type be already normalized.  This option is provided to allow
 calling this function repeatedly without normalizing the same weights
 vector each time.

 The following are more *<optional keyword arguments*>:
 - `init` is a matrix to be used as initialization for the mean. If no matrix is provided, the instance of [generalized means](@ref) with parameter ``p`` will be used.
 - `tol` is the tolerance for the convergence. The smaller this number (it must be positive) the closer the algorithm gets to the true solution;
 - if `⍰` is true, the convergence attained at each iteration is printed.

!!! note "Nota Bene"
  In normal circumstances this algorithm converges monothonically.
  If the algorithm diverges a **warning** is printed indicating the iteration
  when this happened and the algorithm is interrupted.

 **See**: [modified Bhattacharyya mean](@ref)

 ## Examples
    using LinearAlgebra, PosDefManifold
    # Generate a set of 4 random 3x3 SPD matrices
    ℘=randP(3, 4)

    # unweighted mean
    G, iter, conv = powerMean(℘, 0.5)

    # weights vector, does not need to be normalized
    weights=[1, 2, 3, 1]

    # weighted mean
    G, iter, conv = powerMean(℘, 0.5, w=weights)

    # print the convergence at all iterations
    G, iter, conv = powerMean(℘, 0.5, w=weights, ⍰=true)

    # now suppose ℘ has changed a bit, initialize with G to hasten convergence
    ℘[1]=ℍ(℘[1]+(randP(3)/100))
    G, iter, conv = powerMean(℘, 0.5, w=weights, ⍰=true, init=G)

"""
function powerMean(℘, p::Real;     w::Vector=[], ✓w::Bool=true, init=nothing,
                                    tol=1e-9, ⍰=false)
  if !(-1<=p<=1) @error("The parameter p for power means must be in range [-1...1]")
  else
    if     p ≈-1
            return (meanP(℘, InvEuclidean, w=w, ✓w=✓w), 1, 0)
    elseif p ≈ 0
            LE=meanP(℘, logEuclidean, w=w, ✓w=✓w)
            P, iter1, conv1=powerMean(℘,  0.01, w=w, ✓w=✓w, init=LE, tol=tol, ⍰=⍰)
            Q, iter2, conv2=powerMean(℘, -0.01, w=w, ✓w=✓w, init=P, tol=tol, ⍰=⍰)
            return (geodesic(P, Q,  0.5,  Fisher), iter1+iter2, (conv1+conv2)/2)
    elseif p ≈ 1
                return (meanP(℘, Euclidean, w=w, ✓w=✓w), 1, 0)
    else
        # Set Parameters
        maxIter=500
        n, k = Attributes(℘)
        sqrtn=√n
        absp=abs(p)
        r=-0.375/absp
        w≠[] ? v = GetWeights(w, ✓w, k) : v=[]
        init == nothing ? M = generalizedMean(℘, p; w=v, ✓w=false) : M = ℍ(init)
        p<0 ? X=ℍ(M^(0.5)) : X=ℍ(M^(-0.5))
        X◇, H = similar(X, eltype(X))
        𝒫=similar(℘, eltype(℘))
        if p<0 𝒫=[ℍ(inv(P)) for P in ℘] else 𝒫=℘ end
        iter = 1
        conv = 0.; oldconv=maxpos
        ⍰ && @info("Iterating powerMean Fixed-Point...")

        @inbounds while true
            if isempty(w)
                H=ℍ(𝛍(pow(ℍ(X*P*X), absp) for P in 𝒫))
            else
                H=ℍ(𝚺(ω * pow(ℍ(X*P*X), absp) for (ω, P) in zip(v, 𝒫)))
            end
            X◇=(pow(H, r))*X
            conv=norm(H-I)/sqrtn # relative difference to identity
            ⍰ &&  println("iteration: ", iter, "; convergence: ", conv)
            diverging = conv > oldconv
            diverging ? @warn("powerMean diverged at:", iter) : oldconv=conv
            iter==maxIter || diverging || conv <= tol ? break : X = X◇
            iter += 1
        end # while
    end # if

    if p<0  return ( ℍ((X◇)'*X◇), iter, conv )
    else    return ( ℍ(inv((X◇)'*X◇)), iter, conv ) end
  end # if !(-1<=p<=1)
end



"""
    (1) meanP(P::ℍ, Q::ℍ, metric::Metric=Fisher)
    (2) meanP(℘, metric::Metric=Fisher; <w::Vector=[], ✓w::Bool=true>)

 (1) Mean of two positive definite matrices, passed in arbitrary order as
 arguments ``P`` and ``Q``, using the specified `metric` of type
 [Metric::Enumerated type](@ref). By defult the [Fisher](@ref) metric is used.
 The order is arbitrary as all metrics implemented in **PosDefManifold** are symmetric.
 This is the midpoint of the geodesic.
 For the weighted mean of two positive definite matrices use instead
 the [`geodesic`](@ref) function.
 ``P`` and ``Q`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 (2) [Fréchet mean](@ref) of an 1d array ``℘`` of ``k`` positive definite matrices``{P_1,...,P_k}``,
 with optional non-negative real weights ``w={w_1,...,w_k}`` using the specified
 `metric`as in (1).

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

 **See**: [geodesic](@ref), [mean](@ref).

 ## Examples
    using LinearAlgebra, Statistics, PosDefManifold
    # Generate 2 random 3x3 SPD matrices
    P=randP(3)
    Q=randP(3)
    M=meanP(P, Q, logdet0) # (1)
    M=meanP(P, Q) # (1), uses Fisher metric

    # Generate a set of 4 random 3x3 SPD matrices
    ℘=randP(3, 4)
    # weights vector, does not need to be normalized
    weights=[1, 2, 3, 1]
    M=meanP(℘, Euclidean, w=weights) # (2) weighted Euclidean mean
    M=meanP(℘, Wasserstein)  # (2) unweighted Wassertein mean

"""
meanP(P::ℍ, Q::ℍ, metric::Metric=Fisher) = geodesic(P, Q, 0.5, metric)

function meanP(℘, metric::Metric=Fisher;    w::Vector=[], ✓w::Bool=true)
    # iterative solutions
    if      metric == Fisher
            (G, iter, conv)=powerMean(℘, 0; w=w, ✓w=✓w)
            return G
    elseif  metric == logdet0
            (G, iter, conv)=logdet0Mean(℘; w=w, ✓w=✓w)
            return G
    elseif  metric == Wasserstein
            (G, iter, conv)=wasMean(℘; w=w, ✓w=✓w)
            return G
    end

    # closed-form expressions
    n, k = Attributes(℘)
    isempty(w) ? DoNothing : v = GetWeights(w, ✓w, k)
    if  metric == Euclidean
        if isempty(w)   return ℍ(𝛍(P for P in ℘))
        else            return ℍ(𝚺(ω*P for (ω, P) in zip(v, ℘)))
        end

    elseif metric == invEuclidean
        if isempty(w)   return ℍ(inv(ℍ(𝛍(inv(P) for P in ℘))))
        else            return ℍ(inv(ℍ(𝚺(ω*inv(P) for (ω, P) in zip(v, ℘)))))
        end

    elseif metric == logEuclidean
        if isempty(w)   return ℍ(exp(ℍ(𝛍(log(P) for P in ℘))))
        else            return ℍ(exp(ℍ(𝚺(ω*log(P) for (ω, P) in zip(v, ℘)))))
        end

    elseif metric == ChoEuclidean
        if isempty(w)   L = 𝛍(choL(P) for P in ℘)
        else            L = 𝚺(ω*choL(P) for (ω, P) in zip(v, ℘))
        end
        return ℍ(L*L')

    elseif metric == logCholesky # Aggiusta!
        L℘=[choL(P) for P in ℘]
        #if ω==0 L=mean(tril(L℘[i],-1)      for i in 1:k) + exp(mean(𝑓𝑫(L℘[i], log)      for i in 1:k))
        if isempty(w)
            L=𝛍(tril(L℘[i],-1) for i in 1:k)
            for l in 1:n L[l, l] = exp(𝛍(log(L℘[i][l, l]) for i in 1:k)) end
        #else    L=𝚺(ζ[i]*tril(L℘[i],-1) for i in 1:k)/k + exp(𝚺(𝑓𝑫(ζ[i]*L℘[i], log) for i in 1:k))/k end
        else
            L=𝛍(v[i]*tril(L℘[i],-1) for i in 1:k)
            for l in 1:n L[l, l] = exp(𝚺(v[i]*log(L℘[i][l, l]) for i in 1:k)) end
        end
        return ℍ(L*L')

    elseif metric == Jeffrey
        P=meanP(℘, Euclidean; w=w, ✓w=✓w)
        Q=meanP(℘, invEuclidean; w=w, ✓w=✓w)
        P½, P⁻½=pow(P, 0.5, -0.5)
        return ℍ(P½ * sqrt(ℍ(P⁻½ * Q * P⁻½)) * P½)

    elseif metric == VonNeumann
        @warn "function RiemannianGeometryP.meanP and .geodesic not defined for metric $metric"

    else
        @warn "in RiemannianGeometryP.meanP function: the chosen 'metric' does not exist"
    end # if metric
end # function


# -----------------------------------------------------------
# 5. Tangent Space Tools
# -----------------------------------------------------------

"""
    logMap(P::ℍ, G::ℍ, metric::Metric=Fisher)

 *Logaritmic Map:* map a positive definite matrix ``P`` from the SPD or
 Hermitian manifold into the tangent space at base-point ``G`` using the [Fisher](@ref) metric.

 ``P`` and ``G`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 The map is defined as

 `` Log_G(P)=S=G^{1/2}\\textrm{log}\\big(G^{-1/2}PG^{-1/2}\\big)G^{1/2}``.

 The result is an `Hermitian` matrix.
 The inverse operation is [`expMap`](@ref).

 **Arguments** `(P, G, metric)`:
 - ``P`` is the positive definite matrix to be projected onto the tangent space,
 - ``G`` is the tangent space base point,
 - `metric` is a metric of type [Metric::Enumerated type](@ref).

 Currently only the [Fisher](@ref) metric is supported for tangent space operations.

 **See also**: [`vecP`](@ref)

 ## Examples
    using PosDefManifold
    P=randP(3)
    Q=randP(3)
    G=meanP(P, Q)
    # projecting P at the base point given by the geometric mean of P and Q
    S=logMap(P, G)
"""
function logMap(P::ℍ, G::ℍ, metric::Metric=Fisher)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * log(G⁻½ * P * G⁻½') * G½')
    else @warn "in RiemannianGeometryP.logMap function:
                 only the Fisher metric is supported for the logarithmic map."
    end
end

"""

    expMap(S::ℍ, G::ℍ, metric::Metric=Fisher)

 *Exponential Map:* map an `Hermitian` matrix ``S`` from the tangent space at base
 point ``G`` into the SPD or Hermitian manifold (using the [Fisher](@ref) metric).

 ``S`` and ``G`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

 The map is defined as

 `` Exp_G(S)=P=G^{1/2}\\textrm{exp}\\big(G^{-1/2}SG^{-1/2}\\big)G^{1/2}``.

 The result is a positive definite matrix.
 The inverse operation is [`logMap`](@ref).

 **Arguments** `(S, G, metric)`:
 - ``S`` is a Hermitian matrix, real or complex, to be projected on the SPD or Hermitian manifold,
 - ``G`` is the tangent space base point,
 - `metric` is a metric of type [Metric::Enumerated type](@ref).

  Currently only the Fisher metric is supported for tangent space operations.

 ## Examples
    using PosDefManifold, LinearAlgebra
    P=randP(3)
    Q=randP(3)
    G=meanP(P, Q, Fisher)
    # projecting P on the tangent space at the Fisher mean base point G
    S=logMap(P, G)
    # adding the identity in the tangent space and reprojecting back onto the manifold
    H=expMap(ℍ(S+I), G)
"""
function expMap(S::ℍ, G::ℍ, metric::Metric=Fisher)
    if   metric==Fisher
         G½, G⁻½=pow(G, 0.5, -0.5)
         return ℍ(G½ * exp(G⁻½ * S * G⁻½') * G½')
    else @warn "in RiemannianGeometryP.expMap function:
              only the Fisher metric is supported for the exponential map"
    end
end


"""
    vecP(S::ℍ)

 *Vectorize* a tangent vector (matrix) ``S`` (*i.e.*, an `Hermitian` matrix):  mat -> vec.

 It gives weight ``1`` to diagonal elements and √2 to off-diagonal elements
 (Barachant et al., 2012)[🎓](@ref).

 The result is a vector holding ``n(n+1)/2`` elements, where ``n``
 is the size of ``S``.

 ``S`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 The inverse operation is provided by [`matP`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(3)
    Q=randP(3)
    G=meanP(P, Q, Fisher)
    # projecting P at the base point given by the geometric mean of P and Q
    S=logMap(P, G)
    # vectorize S
    ς=vecP(S)
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
    G=meanP(P, Q, Fisher)
    # projecting P at onto the tangent space at the Fisher mean base point
    S=logMap(P, G)
    # vectorize S
    ς=vecP(S)
    # Rotate the vector by an orthogonal matrix
    n=Int(size(S, 1)*(size(S, 1)+1)/2)
    U=randP(n)
    v=U*ς
    # Get the point in the tangent space
    S=matP(v)
"""
function matP(ς::Vector)
  n=Int((-1+√(1+8*length(ς)))/2) # Size of the matrix whose vectorization vector v has size length(v)
  S=Matrix{eltype(ς)}(undef, n, n)
  l=0;
  for j in 1:n-1
    l=l+1
    @inbounds S[j, j]=ς[l]
    for i in j+1:n
      l=l+1
      @inbounds S[i, j]=invsqrt2*ς[l];  S[j, i]=S[i, j]
    end
  end
  S[n, n]=ς[end]
  return ℍ(S)
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
