#   Unit linearAlgebra.jl, part of PosDefManifold Package for julia language
#   v 0.1.3 - last update 28th of April 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DESCRIPTION
#   This Unit implements Linear Algebra methods and procedures
#   useful in relation to the Riemannian Geometry on the manifold
#   of Symmetric Positive Definite (SPD) or Hermitian matrices.
#
#   CONTENT
#   1. Matrix normalizations
#   2. Boolean functions of matrices
#   3. Scalar functions of matrices
#   4. Diagonal functions of matrices
#   5. Unitary functions of matrices
#   6. Matrix function of matrices
#   7. Spectral decompositions of positive matrices
#   8. Decompositions involving triangular matrices
# __________________________________________________________________

#  ------------------------
## 1. Matrix Normalizations
#  ------------------------
"""
    det1!(P::ℍ)

 Given a positive definite matrix ``P``, return the best approximant to
 ``P`` from the set of matrices in the [special linear group](https://bit.ly/2W5jDZ6),
 i.e., the closer matrix having det=1. See Bhatia and Jain (2014)[🎓].

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).
 However a catch-all method is defined.

 **See** [det](https://bit.ly/2Y4MnTF).

 **See also**: [`tr1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=det1(P)
    det(Q) # must be 1

"""
det1(P::ℍ) = ℍ(triu(P)/det(P)^(1/size(P, 1)))

det1(P) = P/det(P)^(1/size(P, 1))


"""
    tr1(P::ℍ)

 Given a positive definite matrix ``P``, return the trace-normalized ``P``
 (trace=1).

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).
 However a catch-all method is defined.

 **See**: [Julia trace function](https://bit.ly/2HoOLiM).

 **See also**: [`tr`](@ref), [`det1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=tr1(P)
    tr(Q)  # must be 1

"""
tr1(P::ℍ) = ℍ(triu(P)/tr(P))

tr1(P) = P/tr(P)


"""
    (1) normalizeCol!(X::Matrix, j::Int)
    (2) normalizeCol!(X::Matrix, j::Int, by::Number)
    (3) normalizeCol!(X::Matrix, range::UnitRange)
    (4) normalizeCol!(X::Matrix, range::UnitRange, by::Number)

 Given a general matrix ``X``,
 - (1) normalize the ``j^{th}``column
 - (2) divide the elements of the ``j^{th}`` column by number ``by``
 - (3) normalize the columns in ``range``
 - (4) divide the elements of columns in ``range``  by number ``by``.

 ``by`` is a number of abstract supertype [Number](https://bit.ly/2JwXjGr).
 It should be an integer, real or complex number.

 ``range`` is a [UnitRange](https://bit.ly/2HSfK5J) type.

  No range check nor type check is performed.
  A catch-all method is defined, but keep in mind that Julia
  does not allow normalizing the columns of `Hermitian` matrices.
  (see [typecasting matrices](@ref)).

  **See** [norm](https://bit.ly/2TaAkR0) and [randn](https://bit.ly/2I1Vgrg) for the example

  **See also**: [`colNorm`](@ref), [`colProd`](@ref).

  ## Examples
    using PosDefManifold
    X=randn(10, 20)
    normalizeCol!(X, 2)                  # (1) normalize columns 2
    normalizeCol!(X, 2, 10.0)            # (2) divide columns 2 by 10.0
    normalizeCol!(X, 2:4)                # (3) normalize columns 2 to 4
    X=randn(ComplexF64, 10, 20)
    normalizeCol!(X, 3)                  # (1) normalize columns 3
    normalizeCol!(X, 3:6, (2.0 + 0.5im)) # (4) divide columns 3 to 5 by (2.0 + 0.5im)

"""
function normalizeCol!(X::Matrix{T}, j::Int) where T<:RealOrComplex
    w=colNorm(X, j)
    for i=1:size(X, 1) @inbounds X[i, j]/=w end
end

function normalizeCol!(X, j::Int)
    w=colNorm(X, j)
    for i=1:size(X, 1) @inbounds X[i, j]/=w end
end

normalizeCol!(X::Matrix{T}, j::Int, by::Number) where T<:RealOrComplex = for i=1:size(X, 1) @inbounds X[i, j]/=by end

normalizeCol!(X, j::Int, by::Number) = for i=1:size(X, 1) @inbounds X[i, j]/=by end

normalizeCol!(X::Matrix{T}, range::UnitRange) where T<:RealOrComplex = for j in range normalizeCol!(X, j) end

normalizeCol!(X, range::UnitRange) = for j in range normalizeCol!(X, j) end

normalizeCol!(X::Matrix{T}, range::UnitRange, by::Number) where T<:RealOrComplex = for j in range normalizeCol!(X, j, by) end

normalizeCol!(X, range::UnitRange, by::Number) = for j in range normalizeCol!(X, j, by) end


#  -------------------------------
## 2. Boolean Functions of Matrices
#  -------------------------------

"""
```
    (1) ispos(   λ::Vector{T};
                <tol::Real=0, rev=true, 🔔=true, msg="">) where T<:Real

    (2) ispos(   Λ::Diagonal{T};
                <tol::Real=0, rev=true, 🔔=true, msg="">) where T<:Real
```

 Return ``true`` if all numbers in (1) real vector ``λ`` or in (2) real diagonal
 matrix ``Λ`` are not inferior to ``tol``, otherwise return ``false``. This may be used,
 for example, in spectral functions to check that all eigenvalues are positive.
 ``tol`` defaults to the square root of `Base.eps` of the type of ``λ`` (1)
 or ``Λ`` (2). This corresponds to requiring positivity beyond about half of
 the significant digits.

 The following are *<optional keyword arguments>*:

 - If ``rev=true`` the (1) elements in ``λ`` or (2) the diagonal elements in ``Λ`` will be chacked in reverse order.
 This is done for allowing a very fast
 check when the elements are sorted where to start checking.

 If the result is ``false``:
 - if ``🔔=true`` a bell character will be printed. In most systems this will ring a bell on the computer.
 - if string ``msg`` is provided, a warning will print ``msg`` followed by:
 "at position *pos*", where *pos* is the position where the
 first non-positive element has been found.

```
 ## Examples
 using PosDefManifold
 a=[1, 0, 2, 8]
 ispos(a, msg="non-positive element found")

 # it will print:
 # ┌ Warning: non-positive element found at position 2
 # └ @ [here julie will point to the line of code issuing the warning]
```
 """
function ispos( λ::Vector{T};   tol::Real=0, rev=true, 🔔=true, msg="")  where T<:Real
    tol==0 ? tolerance = √eps(eltype(λ)) : tolerance = tol
    rev ? iterations = (length(λ):-1:1) : iterations=(1:length(λ))
    for i in iterations
        if λ[i]<tolerance
            🔔 && print('\a') # print('\a') sounds a bell
            length(msg)>0 && @warn("function ispos(linearAlgebra.jl) "*msg* " at position $i")
            return false; break
        end
    end
    return true
end

function ispos( Λ::Diagonal{T};   tol::Real=0, rev=true, 🔔=true, msg="")  where T<:Real
    return ispos( diag(Λ); tol=tol, rev=rev, 🔔=🔔, msg=msg)
end


#  -------------------------------
## 3. Scalar functions of matrices
#  -------------------------------

"""
    (1) colProd(X::Matrix, j::Int, l::Int)
    (2) colProd(X::Matrix, j::Int, l::Int)

 (1) Given a general matrix ``X``, comprised of real or complex elements,
 return the dot product of the ``j^{th}`` and ``l^{th}`` columns, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*x_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and ``^*`` the complex conjugate.

 (2) Given two general matrices ``X`` and ``Y``, comprised of real or complex elements,
 return the dot product of the ``j^{th}`` column of ``X`` and the ``l^{th}`` column
 of ``Y``, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*y_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and of ``Y`` and ``^*`` the complex conjugate.

 ``X`` and of ``Y`` may have a different number of columns.
 A catch-all method is defined.

 Arguments ``j`` and ``l`` must be positive integers in range
 (1) `j,l in 1:size(X, 2)` and (2) `j in 1:size(X, 2), l in 1:size(Y, 2)`.

 **See also**: [`normalizeCol!`](@ref), [`colNorm`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    p=colProd(X, 1, 3)
    Y=randn(10, 30)
    q=colProd(X, Y, 2, 25)

"""
colProd(X::Matrix{T}, j::Int, l::Int) where T<:RealOrComplex =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], X[:, l]))

colProd(X, j::Int, l::Int) =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], X[:, l]))

colProd(X::Matrix{T}, Y::Matrix{T}, j::Int, l::Int) where T<:RealOrComplex =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], Y[:, l]))

colProd(X, Y, j::Int, l::Int) =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], Y[:, l]))


"""
    colNorm(X::Matrix, j::Int)

 Return the Euclidean norm of the ``j^{th}`` column of general matrix ``X``.
 No range check nor type check is performed. A catch-all method is defined.

 **See also**: [`normalizeCol!`](@ref), [`colProd`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    normOfSecondColumn=colNorm(X, 2)

"""
colNorm(X::Matrix{T}, j::Int) where T<:RealOrComplex = √sumOfSqr(X, j)

colNorm(X, j::Int) = √sumOfSqr(X, j)


"""
    (1) sumOfSqr(A::Array)
    (2) sumOfSqr(X::Matrix, j::Int)
    (3) sumOfSqr(X::Matrix, range::UnitRange)

 Return
 - (1) the sum of square of the elements in an array ``A`` of any dimensions.
 - (2) the sum of square of the ``j^{th}`` column of a matrix ``X``.
 - (3) the sum of square of the columns of ``X`` in a given range.

 Note that only (1) works for arrays of any dimensions and that
 if ``A`` is a matrix (1) returns the square of the [Frobenius norm](https://bit.ly/2Fi10eH):
 ``\\sum |a_{ij}|^2. ``

 For (1), if ``A`` is a matrix flagged by Julia as `Hermtian` or as
 `LowerTriangular`, only the lower triangular part of ``A`` is used.

**Arguments**

 (1)  `(A)`:
 - ``A`` is an array of any dimensions (e.g., a vector, matrix or tensor), real or complex.

 (2) `(X, j)`:
 - ``X`` is a generic matrix, real or complex;
 - ``j`` is a positive integer in range `1:size(X, 2)`.

 (3) `(X, range)`:
 - ``X`` is a generic matrix, real or complex;;
 - ``range`` is a [UnitRange type](https://bit.ly/2HDoFbk).

 **See also**: [`sumOfSqrDiag`](@ref), [`sumOfSqrTril`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    sum2=sumOfSqr(X)        # (1) sum of squares of all elements
    sum2=sumOfSqr(X, 1)     # (2) sum of squares of elements in column 1
    sum2=sumOfSqr(X, 2:4)   # (3) sum of squares of elements in column 2 to 4

"""
sumOfSqr(A::Array{T}) where T<:RealOrComplex = 𝚺(abs2(a) for a in A)

function sumOfSqr(H::Union{ℍ, 𝕃})
    r=size(H, 1)
    s=eltype(H)(0)
    for j=1:size(H, 2)-1
        @inbounds s+=abs2(H[j, j])
        for i=j+1:r @inbounds s+=2*abs2(H[i, j]) end
    end
    @inbounds s+=abs2(H[r, r])
    return s
end

sumOfSqr(X::Matrix{T}, j::Int) where T<:RealOrComplex = 𝚺(abs2.(X[:, j]))

sumOfSqr(H::ℍ, j::Int) = 𝚺(abs2.(H[:, j]))

sumOfSqr(X::Matrix{T}, range::UnitRange) where T<:RealOrComplex =
         𝚺(sumOfSqr(X, j) for j in range)

sumOfSqr(H::ℍ, range::UnitRange) = 𝚺(sumOfSqr(H, j) for j in range)


"""
    (1) sumOfSqrDiag(X::Matrix)
    (2) sumOfSqrDiag(D::Diagonal)

 **alias**: `ssd`

 Return (1) the sum of squares of the diagonal elements in general matrix ``X``
 comprised of real or complex numbers.
 If ``X`` is rectangular, the main diagonal is considered.

 It also return (2) the sum of squares of real diagonal matrix ``Λ``.

 No range check nor type check is performed. A catch-all method is defined.

 **See also**: [`sumOfSqr`](@ref), [`sumOfSqrTril`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    X=randn(10, 20)
    sumDiag2=sumOfSqrDiag(X) # (1)
    sumDiag2=sumOfSqrDiag(𝔻(X)) # (2) 𝔻=LinearAlgebra.Diagonal

"""
sumOfSqrDiag(X::Matrix{T}) where T<:RealOrComplex =
             𝚺(abs2(X[i, i]) for i=1:minimum(size(X)))

sumOfSqrDiag(Λ::Diagonal) = 𝚺(abs2(Λ[i, i]) for i=1:size(Λ, 1))

sumOfSqrDiag(X) = 𝚺(abs2(X[i, i]) for i=1:minimum(size(X)))

ssd=sumOfSqrDiag

"""
    sumOfSqrTril(X::Matrix, k::Int=0)

**alias**: `sst`

 Given a general matrix ``X``, return the sum of squares of the elements
 in the lower triangle ``X`` up to the ``k^{th}`` underdiagonal.

 ``X`` may be rectangular. ``k`` must be in range `1-size(X, 1):0`.

  See julia [tril(M, k::Integer)](https://bit.ly/2Tbx8o7) function
 for numbering of diagonals.
  No range check nor type check is performed. A catch-all method is defined.

 **See also**: [`sumOfSqr`](@ref), [`sumOfSqrDiag`](@ref).

 ## Examples
    using PosDefManifold
    A=[4. 3.; 2. 5.; 1. 2.]
    #3×2 Array{Float64,2}:
    # 4.0  3.0
    # 2.0  5.0
    # 1.0  2.0

    s=sumOfSqrTril(A, -1)
    # 9.0 = 1²+2²+2²

    s=sumOfSqrTril(A, 0)
    # 50.0 = 1²+2²+2²+4²+5²

"""
function sumOfSqrTril(X::Matrix{T}, k::Int=0) where T<:RealOrComplex
    (r, c)=size(X)
    if k<(1-r) || k>(c-1)
        @warn "in LinearAmgebraInP.sumOfSqrTRil function: argument k is out of bounds"
    else
        s=0.0; @inbounds for j=1:c, i=max(j-k, 1):r s+=abs2(X[i, j]) end
        return s
    end
end

function sumOfSqrTril(X, k::Int=0)
    (r, c)=size(X)
    if k<(1-r) || k>(c-1)
        @warn "in LinearAmgebraInP.sumOfSqrTRil function (catch-all input): argument k is out of bounds"
    else
        s=0.0; @inbounds for j=1:c, i=max(j-k, 1):r s+=abs2(X[i, j]) end
        return s
    end
end

sst=sumOfSqrTril

"""
    (1) tr(P::ℍ, Q::ℍ)
    (2) tr(P::ℍ, Q::Matrix)

 Given (1) two positive definite matrix ``P`` and ``Q``,
 return the trace of the product ``PQ``.
 This is real even if ``P`` and ``Q`` are complex.

 ``P`` must always be flagged as `Hermitian`. See [typecasting matrices](@ref).
 In (2) ``Q`` is a generic `Matrix` object,
 in which case return
 - a real trace if the product ``PQ`` is real or if it has all positive real eigenvalues.
 - a complex trace if the product ``PQ`` is not real and has complex eigenvalues.

 ## Math
 Let ``P`` and ``Q`` be `Hermitian` matrices, using the properties of the trace
 (e.g., the cyclic property and the similarity invariance) you can use this
 function to fast compute the trace of many expressions. For example:

 ``\\textrm{tr}(PQ)=\\textrm{tr}(P^{1/2}QP^{1/2})``

 and

 ``\\textrm{tr}(PQP)=\\textrm{tr}(P^{2}Q)`` (see example below).


 **See**: [trace](https://bit.ly/2HoOLiM).

 **See also**: [`tr1`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
    Q=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
    tr(P, Q) ≈ tr(P*Q) ? println(" ⭐ ") : println(" ⛔ ")
    tr(P, Q) ≈ tr(sqrt(P)*Q*sqrt(P)) ? println(" ⭐ ") : println(" ⛔ ")
    tr(sqr(P), Q) ≈ tr(P*Q*P) ? println(" ⭐ ") : println(" ⛔ ")

"""
function tr(P::ℍ, Q::ℍ)
    a = 𝚺(colProd(P, Q, i, i) for i=1:size(P, 2))
    if real(a)<0 return a else return real(a) end;
end

function tr(P::ℍ, Q::Matrix)
    λ = [colProd(P, Q, i, i) for i=1:size(P, 2)]
    OK=true
    for l in λ
        if imag(l) ≉  0
            OK=false
            break
        end
    end
    if OK return real(𝚺(λ)) else return 𝚺(λ) end
end


"""
    quadraticForm(v::Vector{T}, X::Matrix{T}) where T<:RealOrComplex

 **alias**: `qufo`

 Given a vector ``v`` and a matrix ``X``, compute the quadratic form

 ``v^*Xv``.

 ``v`` and ``X`` may be real or complex. If they are real, only the
 lower triangular part of ``X`` is used.

 ``X`` may be a generic `Matrix`, may be flagged by Julia as an `Hermitian`
 matrix or as a `LowerTriangular` matrix. In the latter case it must be real.

 ## Math

 For ``v`` and ``X`` real, the quadratic form is

 ```sum_i(v_i^2x_{ii})+sum_i>j(2v_iv_jx_{ij})``.


 ## Examples
    using PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    v=randn(5)
    q1=quadraticForm(v, P) # or q1=quad(v, P)
    # obtain a lower Triangular view of P
    L=LowerTriangular(Matrix(P)) # or L=𝕃(Matrix(P))
    q2=quadraticForm(v, L)
    q1 ≈ q2 ? println(" ⭐ ") : println(" ⛔ ")

"""
function quadraticForm(v::Vector{T}, X::Matrix{T}) where T<:RealOrComplex
    if T<:Real
        r=length(v)
        s=0.
        for i=1:r @inbounds s+=(v[i]^2 * X[i, i]) end # Diagonal
        for j=1:r-1, i=j+1:r @inbounds s+=2*v[i]*v[j]*X[i, j] end # Off-diagonal
        return s
    else
        return v'*X*v
    end
end

quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:RealOrComplex =
              quadraticForm(v, Matrix(P))

quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real =
              quadraticForm(v, Matrix(L))

qufo=quadraticForm


"""
    fidelity(P::ℍ, Q::ℍ)

 Given two positive definte matrices ``P`` and ``Q``, return their *fidelity*:

  ``tr\\big(P^{1/2}QP^{1/2}\\big)^{1/2}.``

  This is used in quantum physics and is related to the
 [Wasserstein](@ref) metric. See for example Bhatia, Jain and Lim (2019b)[🎓](@ref).

 ``P`` and ``Q`` must be flagged as `Hermitian`.
 See [typecasting matrices](@ref),
 however a catch-all method is defined.

 ## Examples
    using PosDefManifold
    P=randP(5);
    Q=randP(5);
    f=fidelity(P, Q)

"""
function fidelity(P::ℍ, Q::ℍ)
    A = √(P)
    return tr(√ℍ(A*Q*A'))
end


#  ---------------------------------
## 4. Diagonal functions of matrices
#  ---------------------------------
"""
    (1) fDiagonal(func::Function, X::ℍ, k::Int=0)
    (2) fDiagonal(func::Function, X::𝕃, k::Int=0)
    (3) fDiagonal(func::Function, X::Diagonal, k::Int=0)
    (4) fDiagonal(func::Function, X::Matrix, k::Int=0)

 **alias**: `𝑓𝔻`

 Applies function `func` element-wise to the elements of the ``k^{th}``
 diagonal of matrix ``X`` (square in (1-3) and of dimension *r⋅c* in (4))
 and return a diagonal matrix with these elements.

 See julia [tril(M, k::Integer)](https://bit.ly/2Tbx8o7) function
 for numbering of diagonals.

 If the matrix is Diagonal (3) `k` must be zero.
 If the matrix is lower triangular (2) `k` cannot be positive.

 Note that if ``X`` is rectangular the dimension of the result depends
 on the size of ``X`` and on the chosen diagonal.
 For example,
 - *r ≠ c* and ``k``=0 (main diagonal), the result will be of dimension min*(r,c)*⋅*min(r,c)*,
 - ``X`` *3⋅4* and ``k=-1``, the result will be *2⋅2*,
 - ``X`` *3⋅4* and ``k=1``, the result will be *3⋅3*, etc.

!!! note "Nota Bene"
    The function `func` must support the `func.` syntax and therefore
    must be able to apply element-wise to the elements of the chosen diagonal
    (this includes anonymous functions). If the input matrix is complex, the function `func`
    must be able to support complex arguments.

 ## Examples
    using PosDefManifold
    P=randP(5) # use P=randP(ComplexF64, 5) for generating an Hermitian matrix
    D=fDiagonal(inv, P, -1) # diagonal matrix with the inverse of the first sub-diagonal of P
    (Λ, U) = evd(P)         # Λ holds the eigenvalues of P, see evd
    Δ=fDiagonal(log, Λ)     # diagonal matrix with the log of the eigenvalues
    Δ=fDiagonal(x->x^2, Λ)  # using an anonymous function for the square of the eigenvalues
"""
fDiagonal(func::Function, D::𝔻, k::Int=0) = 𝔻(func.(D))

function fDiagonal(func::Function, L::𝕃, k::Int=0)
 if k>0 @error("in function fDiagonal (linearAlgebra.jl): k argument cannot be positive.")
 else return 𝔻(func.(diag(L, k)))
 end
end

fDiagonal(func::Function, P::ℍ, k::Int=0) = 𝔻(func.(diag(P, k)))

fDiagonal(func::Function, X::Matrix{T}, k::Int=0) where T<:RealOrComplex= 𝔻(func.(diag(X, k)))

𝑓𝔻=fDiagonal



#  -------------------------------
## 5. Unitary functions of matrices
#  -------------------------------
"""
    mgs(T::Matrix, numCol::Int=0)

 Modified (stabilized) [Gram-Schmidt orthogonalization](https://bit.ly/2YE6zvy)
 of the columns of square or tall matrix ``T``, which can be comprised of real
 or complex elements.
 The orthogonalized ``T`` is returned by the function. ``T`` is not changed.

 All columns are orthogonalized by default. If instead argument `numCol` is provided,
 then only the first `numCol` columns of ``T`` are orthogonalized.
 In this case only the firt `numCol` columns will be returned.

 ## Examples
    using LinearAlgebra, PosDefManifold
    X=randn(10, 10);
    U=mgs(X)        # result is 10⋅10
    U=mgs(X, 3)     # result is 10⋅3
    U'*U ≈ I ? println(" ⭐ ") : println(" ⛔ ")
    # julia undertands also:
    U'U ≈ I ? println(" ⭐ ") : println(" ⛔ ")

"""
function mgs(T::Matrix, numCol::Int=0)
    r=size(T, 1)
    if numCol != 0 && numCol in 2:size(T, 2) c=numCol else c=size(T, 2) end
    U=Matrix{eltype(T)}(undef, r, c)
    U[:, 1] = T[:, 1]/colNorm(T, 1);
    for i in 2:c
        U[:, i]=T[:, i]
        for j in 1:i-1
            s = colProd(U, j, i) / sumOfSqr(U, j)
            U[:, i] -= s * U[:, j]
        end
        normalizeCol!(U, i)
    end
    return U
end # mgs function

#  ------------------------------
## 6. Matrix function of matrices
#  ------------------------------

#  -----------------------------------------------
## 7. Spectral decompositions of positive matrices
#  -----------------------------------------------

"""
    evd(S::ℍ)

 Given a positive semi-definite matrix ``S``,
 returns a 2-tuple ``(Λ, U)``, where ``U`` is the matrix holding in columns
 the eigenvectors and ``Λ`` is the matrix holding the eigenvalues on the diagonal.
 This is the output of Julia `eigen` function in ``UΛU'=S`` form.

 As for the `eigen` function, the eigenvalues and associated
 eigenvectors are sorted by increasing values of eigenvalues.

 ``S`` must be flagged by Julia as `Hermitian`.
 See [typecasting matrices](@ref).

 **See also**: [`spectralFunctions`](@ref).

 ## Examples
    using PosDefManifold
    A=randn(3, 3);
    S=ℍ(A+A');
    Λ, U=evd(S); # which is equivalent to (Λ, U)=evd(P)
    (U*Λ*U') ≈ S ? println(" ⭐ ") : println(" ⛔ ")
    # => UΛU'=S, UΛ=SU, ΛU'=U'S
"""
function evd(S::ℍ) # return tuple (Λ, U)
    F = eigen(S)
    return  𝔻(F.values), F.vectors # 𝔻=LinearAlgebra.Diagonal
end


"""
    spectralFunctions(P::ℍ, func)

 This is the *mother function* for all spectral functions of eigenvalues implemented
 in this library, which are:
 - `pow`     (power),
 - `isqrt`   (inverse square root).

 The function `sqr` (square) does not use it, as it can be obtained more
 efficiently by simple multiplication.

 You can use this function if you need another spectral function of eigenvalues
 besides those and those already implemented in the standard package `LinearAlgebra`.
 In general, you won't call it directly.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 The definition of spectral functions for a positive definite matrix ``P``
 is at it follows:

 ``f\\big(P\\big)=Uf\\big(Λ\\big)U',``

 where ``U`` is the matrix holding in columns the eigenvectors of ``P``,
 ``Λ`` is the matrix holding on diagonal its eigenvalues and ``f`` is
 a function applying element-wise to the eigenvalues.


 **Arguments** `(P, func)`;
 - ``P`` is a positive matrix.
 - `func is the function that will be applied on the eigenvalues

!!! note "Nota Bene"
    The function `func` must support the `func.` syntax and therefore
    must be able to apply element-wise to the eigenvalues
    (those include anonymous functions).

 **See also**: [`evd`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    n=5
    P=randP(n) # P=randP(ComplexF64, 5) to generate an Hermitian complex matrix
    noise=0.1;
    Q=spectralFunctions(P, x->x+noise) # add white noise to the eigenvalues
    tr(Q)-tr(P) ≈ noise*n ? println(" ⭐ ") : println(" ⛔ ")

"""
function spectralFunctions(P::ℍ, func::Function)
    F = eigen(P)
    ispos(F.values, msg="function spectralFunctions: at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return ℍ(F.vectors * 𝔻(func.(F.values)) * F.vectors')
end


"""
    pow(P::ℍ, p)        # one argument
    pow(P::ℍ, args...)  # several arguments

 Given a positive definite matrix ``P``, return the power
 ``P^{r_1}, P^{r_2},...``
 for any number of exponents ``r_1, r_2,...``.
 It returns a tuple of as many elements as arguments passed after ``P``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 **Arguments** `(P, arg1, arg2,...)`
 - ``P`` is a positive matrix.
 - ``arg1, arg2,...`` are real numbers.

 **See also**: [`invsqrt`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5);     # use P=randP(ComplexF64, 5) for generating an Hermitian matrix
    Q=pow(P, 0.5);            # =>  QQ=P
    Q, W=pow(P, 0.5, -0.5);
    W*P*W ≈ I ? println(" ⭐ ") : println(" ⛔ ")
    Q*Q ≈ P ? println(" ⭐ ") : println(" ⛔ ")
    R, S=pow(P, 0.3, 0.7);
    R*S ≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
pow(P::ℍ, p)=spectralFunctions(P, x->x^p) # one argument
function pow(P::ℍ, args...)               # several arguments
    (Λ, U) = evd(P)
    ispos(Λ, msg="function Rpow: at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return  (ℍ(U * Λ^p * U') for p in args)
end



"""
    invsqrt(P::ℍ)

 Given a positive definite matrix ``P``, compute the inverse of the principal
 square root ``P^{-1/2}``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 **See**: [typecasting matrices](@ref).

 **See also**: [`pow`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(ComplexF64, 5);
    Q=invsqrt(P);
    Q*P*Q ≈ I ? println(" ⭐ ") : println(" ⛔ ")

"""
invsqrt(P::ℍ) = spectralFunctions(P, x->1/sqrt(x));


"""
    sqr(P::ℍ)

 Given a positive definite matrix ``P``, compute its square ``P^{2}``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 **See also**: [`pow`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(5);
    P²=sqr(P);  # =>  P²=PP
    sqrt(P²)≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
sqr(P::ℍ) = ℍ(P*P)


"""
    powerIterations(H, q;
            <evalues=false, tol::Real=0, maxiter=300, ⍰=false>)

 **alias**: `powIter`

 Compute the ``q`` eigenvectors associated to the ``q`` largest (real) eigenvalues
 of matrix ``H`` using the [power iterations](https://bit.ly/2JSo0pb) +
 [Gram-Schmidt orthogonalization](https://bit.ly/2YE6zvy) as suggested by Strang.
 The eigenvectors are returned with the same type as ``H``.

 ``H`` must have real eigenvalues. It may be a generic `Matrix` object or may be
 flagged by Julia as `Hermitian`. See [typecasting matrices](@ref).
 In both cases it must be a symmetric matrix if it is real
 or an Hermitian matrix if it is complex. ``H`` may also be `LowerTriangular`,
 but only if it is real (see below).

 The following are *<optional keyword arguments>*:
 - ``tol`` is the tolerance for the convergence of the power method (see below).
 - ``maxiter`` is the maximum number of iterations allowed for the power method
 - if ``⍰=true``, the convergence of all iterations will be printed.
 - if ``evalues=true``, return the 4-tuple ``(Λ, U, iterations, covergence)``
 - if ``evalues=false`` return the 3-tuple ``(U, iterations, covergence)``


!!! note "Nota Bene"
    Differently from the [`evd`](@ref) function, the eigenvectors and
    eigenvalues are sorted by decreasing order of eigenvalues.

    If ``H`` is real, only its lower triangular part is used.
    In this case a BLAS routine is used for computing the power iterations.
    See [BLAS routines](@ref).

    ``tol`` defaults to the square root of `Base.eps` of the type of ``H``.
    This corresponds to requiring equality for the convergence criterion
    over two successive iterations of about half the significant digits.


**See also**: [`mgs`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    # Generate an Hermitian (complex) matrix
    H=randP(ComplexF64, 10);
    # all eigenvectors
    U, iterations, convergence=powIter(H, size(H, 2), ⍰=true)
    # 3 eigenvectors and eigenvalues
    Λ, U, iterations, convergence=powIter(H, 3, evalues=true);
    U'*U ≈ I && U*Λ*U'≈H ? println(" ⭐ ") : println(" ⛔ ")

    # passing a `Matrix` object
    Λ, U, iterations, convergence=powIter(Matrix(H), 3, evalues=true)

    # passing a `LowerTriangular` object (must be a real matrix in this case)
    L=𝕃(randP(10))
    Λ, U, iterations, convergence=powIter(L, 3, evalues=true)

"""
function powerIterations(H::Matrix, q::Int;
                         evalues=false, tol::Real=0, maxiter=300, ⍰=false)

    (n, q², type) = size(H, 1), q^2, eltype(H)
    tol==0 ? tolerance = √eps(real(type)) : tolerance = tol
    msg1="Power Iterations reached a saddle point at:"
    msg2="Power Iterations reached the max number of iterations at:"
    U=randn(type, n, q) # initialization
    normalizeCol!(U, 1:q)
    💡=similar(U) # 💡 is the poweriteration matrix
    (iter, conv, oldconv) = 1, 0., maxpos
    ⍰ && @info("Running Power Iterations...")
    while true
        # power iteration of q vectors and their Gram-Schmidt Orthogonalization
        type<:Real ? 💡=mgs(BLAS.symm('L', 'L', H, U)) : 💡=mgs(H*U)
        conv=norm((💡)' * U-I) / q²
        saddlePoint = conv ≈ oldconv  && @info(msg1, iter)
        overRun     = iter == maxiter && @warn(msg2, iter)
        #diverged    = conv > oldconv && @warn(msg3, iter)
        ⍰ && println("iteration: ", iter, "; convergence: ", conv)
        if conv<=tolerance || saddlePoint==true ||  overRun==true
            break
        else U = 💡 end
        oldconv=conv
        iter += 1
    end # while
    if evalues == false
        return (💡, iter, conv)
    else
        type<:Real ? d=[quadraticForm(💡[:, i], H) for i=1:q] : d=[💡[:, i]'*H*💡[:, i] for i =1:q]
        return (𝔻(real(d)), 💡, iter, conv)
    end
end

powerIterations(H::ℍ, q::Int;
        evalues=false, tol::Real=0, maxiter=300, ⍰=false) =
    powerIterations(Matrix(H), q; evalues=evalues, tol=tol, maxiter=maxiter, ⍰=⍰)

powerIterations(L::𝕃{T}, q::Int;
        evalues=false, tol::Real=0, maxiter=300, ⍰=false) where T<:Real =
    powerIterations(Matrix(L), q; evalues=evalues, tol=tol, maxiter=maxiter, ⍰=⍰)

powIter=powerIterations

#  -----------------------------------------------
## 8. Decompositions involving triangular matrices
#  -----------------------------------------------
"""
    choL(P::Union{ℍ, Matrix})

 Given a positive matrix ``P``, return the *Cholesky lower triangular factor* ``L``
 such that ``LL'=P``. To obtain ``L'`` or both ``L`` and ``L'``, use instead
 julia function [cholesky(P)](https://bit.ly/2u9Hw5P).

 ``P`` sould be flagged as `Hermitian` - see
 [typecasting matrices](@ref) - but a method for generic matrices
 is also provided.

 On output, ``L`` is of type [`LowerTriangular`](https://bit.ly/2U511f3).

 ## Examples
    using PosDefManifold
    P=randP(5);
    L=choL(P);
    L*L'≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
function choL(P::Union{ℍ, Matrix})
    choP = cholesky(P)
    return choP.L
end
