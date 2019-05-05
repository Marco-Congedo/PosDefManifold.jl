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
    det1(X::𝕄)
    det1(X::ℍ)

 Given a real or complex square `Matrix` or `Hermitian` matrix ``X``,
 return the best approximant to
 ``X`` from the set of matrices in the [special linear group](https://bit.ly/2W5jDZ6),
 *i.e.*, the closer matrix having determinant equal to 1.
 See Bhatia and Jain (2014)[🎓].

  **See** [det](https://bit.ly/2Y4MnTF).

 **See also**: [`tr1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=det1(P)
    det(Q) # must be 1

"""
det1(X::ℍ) = ℍ(triu(X)/det(X)^(1/size(X, 1)))

det1(X::𝕄) = X/det(X)^(1/size(X, 1))


"""
    tr1(X::ℍ)
    tr1(X::𝕄)

 Given a real or complex square `Matrix` or `Hermitian` matrix ``X``,
 return the trace-normalized ``X``
 (trace=1).

  **See**: [Julia trace function](https://bit.ly/2HoOLiM).

 **See also**: [`tr`](@ref), [`det1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=tr1(P)
    tr(Q)  # must be 1

"""
tr1(X::ℍ) = ℍ(triu(X)/tr(X))

tr1(X::𝕄) = X/tr(X)


"""
    (1) normalizeCol!(X::𝕄, j::Int)
    (2) normalizeCol!(X::𝕄, j::Int, by::Number)
    (3) normalizeCol!(X::𝕄, range::UnitRange)
    (4) normalizeCol!(X::𝕄, range::UnitRange, by::Number)


 Given a `Matrix` ``X``,
 - (1) normalize the ``j^{th}``column
 - (2) divide the elements of the ``j^{th}`` column by number ``by``
 - (3) normalize the columns in ``range``
 - (4) divide the elements of columns in ``range``  by number ``by``.

 ``by`` is a number of abstract supertype [Number](https://bit.ly/2JwXjGr).
 It should be an integer, real or complex number.

 ``range`` is a [UnitRange](https://bit.ly/2HSfK5J) type.

!!! note "Nota Bene"
    Julia does not allow normalizing the columns of `Hermitian` matrices.
    If you want to call this function for an `Hermitian` matrix see [typecasting matrices](@ref).

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
function normalizeCol!(X::𝕄, j::Int)
    w=colNorm(X, j)
    for i=1:size(X, 1) @inbounds X[i, j]/=w end
end

normalizeCol!(X::𝕄, j::Int, by::Number) =
             for i=1:size(X, 1) @inbounds X[i, j]/=by end

normalizeCol!(X::𝕄, range::UnitRange) =
             for j in range normalizeCol!(X, j) end

normalizeCol!(X::𝕄, range::UnitRange, by::Number) =
             for j in range normalizeCol!(X, j, by) end


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
 matrix ``Λ`` are not inferior to ``tol``, otherwise return ``false``. This is used,
 for example, in spectral functions to check that all eigenvalues are positive.

!!! note "Nota Bene"
    ``tol`` defaults to the square root of `Base.eps` of the type of ``λ`` (1)
     or ``Λ`` (2). This corresponds to requiring positivity beyond about half of
     the significant digits.

 The following are *<optional keyword arguments>*:

 - If ``rev=true`` the (1) elements in ``λ`` or (2) the diagonal elements in ``Λ`` will be chacked in reverse order.
 This is done for allowing a very fast
 check when the elements are sorted and it is known where to start checking.

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
    (1) colProd(X::Union{𝕄, ℍ}, j::Int, l::Int)
    (2) colProd(X::Union{𝕄, ℍ}, Y::Union{𝕄, ℍ}, j::Int, l::Int)

 (1) Given a real or complex `Matrix` or `Hermitian` matrix ``X``,
 return the dot product of the ``j^{th}`` and ``l^{th}`` columns, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*x_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and ``^*`` denotes complex conjugate.

 (2) Given real or complex `Matrix` or `Hermitian` matrices ``X`` and ``Y``,
 return the dot product of the ``j^{th}`` column of ``X`` and the ``l^{th}`` column
 of ``Y``, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*y_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and of ``Y`` and ``^*`` denotes
 the complex conjugate.

!!! note "Nota Bene"
    ``X`` and of ``Y`` may have a different number of columns, but must have
    the same number of rows.

 Arguments ``j`` and ``l`` must be positive integers in range
 - (1) `j,l in 1:size(X, 2)`,
 - (2) `j in 1:size(X, 2), l in 1:size(Y, 2)`.

 **See also**: [`normalizeCol!`](@ref), [`colNorm`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    p=colProd(X, 1, 3)
    Y=randn(10, 30)
    q=colProd(X, Y, 2, 25)

"""
colProd(X::Union{𝕄, ℍ}, j::Int, l::Int) =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], X[:, l]))

colProd(X::Union{𝕄, ℍ}, Y::Union{𝕄, ℍ}, j::Int, l::Int) =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], Y[:, l]))


"""
    colNorm(X::Union{𝕄, ℍ}, j::Int)

 Given a real or complex `Matrix` or `Hermitian` matrix ``X``,
 return the Euclidean norm of its ``j^{th}`` column.

 **See also**: [`normalizeCol!`](@ref), [`colProd`](@ref), [`sumOfSqr`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    normOfSecondColumn=colNorm(X, 2)

"""
colNorm(X::Union{𝕄, ℍ}, j::Int) = √sumOfSqr(X, j)


"""
    (1) sumOfSqr(A::Array)
    (2) sumOfSqr(H::ℍ)
    (3) sumOfSqr(L::𝕃)
    (4) sumOfSqr(D::𝔻)
    (5) sumOfSqr(X::Union{𝕄, ℍ}, j::Int)
    (6) sumOfSqr(X::Union{𝕄, ℍ}, range::UnitRange)

**alias**: `ss`

 Return
 - (1) the sum of square of the elements in an array ``A`` of any dimensions.
 - (2) as in (1), but for an `Hermitian` matrix ``H``, using only the lower triangular part.
 - (3) as in (1), but for a `LowerTriangular` matrix ``L``.
 - (4) as in (1), but for a `Diagonal` matrix ``D`` (sum of squares of diagonal elements).
 - (5) the sum of square of the ``j^{th}`` column of a `Matrix` or `Hermitian` ``X``.
 - (6) the sum of square of the columns of a `Matrix` or `Hermitian` ``X`` in a given range.

 All methods support real and complex matrices.

 Only method (1) works for arrays of any dimensions.

 Methods (1)-(4) return the square of the [Frobenius norm](https://bit.ly/2Fi10eH).

 For method (5), ``j`` is a positive integer in range `1:size(X, 1)`.

 For method (6), ``range`` is a [UnitRange type](https://bit.ly/2HDoFbk).

 **See also**: [`colNorm`](@ref), [`sumOfSqrDiag`](@ref), [`sumOfSqrTril`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    sum2=sumOfSqr(X)        # (1) sum of squares of all elements
    sum2=sumOfSqr(X, 1)     # (2) sum of squares of elements in column 1
    sum2=sumOfSqr(X, 2:4)   # (3) sum of squares of elements in column 2 to 4

"""
sumOfSqr(A::Array) = 𝚺(abs2(a) for a in A)

function sumOfSqr(H::ℍ)
    r=size(H, 1)
    s=real(eltype(H))(0)
    for j=1:size(H, 2)-1
        @inbounds s+=abs2(H[j, j])
        for i=j+1:r @inbounds s+=2*abs2(H[i, j]) end
    end
    @inbounds s+=abs2(H[r, r])
    return s
end

function sumOfSqr(L::𝕃)
    s=real(eltype(L))(0)
    for j=1:size(L, 2), i=j:size(L, 1)
        @inbounds s+=abs2(L[i, j])
    end
    return s
end

sumOfSqr(D::𝔻) = sumOfSqrDiag(D)

sumOfSqr(X::Union{𝕄, ℍ}, j::Int) = 𝚺(abs2.(X[:, j]))

sumOfSqr(X::Union{𝕄, ℍ}, range::UnitRange) = 𝚺(sumOfSqr(X, j) for j in range)

ss=sumOfSqr

"""
    (1) sumOfSqrDiag(X::𝕄)
    (2) sumOfSqrDiag(X::Union{𝔻, ℍ, 𝕃})

 **alias**: `ssd`

 (1) Sum of squares of the diagonal elements in real or complex `Matrix` ``X``.
 If ``X`` is rectangular, the main diagonal is considered.

 (2) Sum of squares of the main diagonal of real or complex `Diagonal`,
 `Hermitian` or `LowerTriangular` matrix ``X``.

 **See also**: [`sumOfSqr`](@ref), [`sumOfSqrTril`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    X=randn(10, 20)
    sumDiag2=sumOfSqrDiag(X) # (1)
    sumDiag2=sumOfSqrDiag(𝔻(X)) # (2) 𝔻=LinearAlgebra.Diagonal

"""
sumOfSqrDiag(X::𝕄) = 𝚺(abs2(X[i, i]) for i=1:minimum(size(X)))

sumOfSqrDiag(X::Union{𝔻, ℍ, 𝕃}) = 𝚺(abs2(X[i, i]) for i=1:size(X, 1))

ssd=sumOfSqrDiag

"""
    sumOfSqrTril(X::Union{𝕄, 𝔻, ℍ, 𝕃}, k::Int=0)

**alias**: `sst`

 Given a real or complex `Matrix`, `Diagonal`, `Hermitian` or
 `LowerTriangular` matrix ``X``,
 return the sum of squares of the elements
 in its lower triangle up to the ``k^{th}`` underdiagonal.

 `Matrix` ``X`` may be rectangular.

 ``k`` must be in range
 - `1-size(X, 1):c-1` for ``X`` `Matrix`, `Diagonal` or `Hermitian`,
 - `1-size(X, 1):0` for ``X`` `LowerTriangular`.

 See julia [tril(M, k::Integer)](https://bit.ly/2Tbx8o7) function
 for numbering of diagonals.

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
function sumOfSqrTril(X::Union{𝕄, 𝔻, ℍ, 𝕃}, k::Int=0)
    (r, c) = size(X)
    X isa(𝕃) ? range = (1-r:0) : range = (1-r:c-1)
    if k in range
        s=eltype(X)(0)
        @inbounds for j=1:c, i=max(j-k, 1):r s+=abs2(X[i, j]) end
        return s
    else
        @error "in LinearAmgebraInP.sumOfSqrTRil function: argument k is out of bounds"
    end
end
sst=sumOfSqrTril

"""
    (1) tr(P::ℍ, Q::ℍ)
    (2) tr(P::ℍ, Q::𝕄)
    (3) tr(D::𝔻{T}, H::Union{ℍ, 𝕄}) where T<:Real
    (4) tr(H::Union{ℍ, 𝕄}, D::𝔻{T}) where T<:Real

 Given (1) two `Hermitian` positive definite matrix ``P`` and ``Q``,
 return the trace of the product ``PQ``.
 This is real even if ``P`` and ``Q`` are complex.

 ``P`` must always be flagged as `Hermitian`. See [typecasting matrices](@ref).
 In (2) ``Q`` is a `Matrix` object,
 in which case return
 - a real trace if the product ``PQ`` is real or if it has all positive real eigenvalues.
 - a complex trace if the product ``PQ`` is not real and has complex eigenvalues.

 Methods (3) and (4) return the trace of the product ``DH`` or ``HD``,
 where ``D`` is a real `Diagonal`matrix and ``H`` an ``Hermitian``
 or ``Matrix`` object. The result if of the same type as the elements of ``H``.

 ## Math
 Let ``P`` and ``Q`` be `Hermitian` matrices, using the properties of the trace
 (e.g., the cyclic property and the similarity invariance) you can use this
 function to fast compute the trace of several expressions. For example:

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
tr(P::ℍ, Q::ℍ) = real(𝚺(colProd(P, Q, i, i) for i=1:size(P, 2)))

function tr(P::ℍ, Q::𝕄)
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


function tr(D::𝔻{T}, H::Union{ℍ, 𝕄}) where T<:Real
    s=eltype(H)(0)
    for i=1:size(D, 1) @inbounds s += D[i, i] * H[i, i] end
    return s
end

tr(H::Union{ℍ, 𝕄}, D::𝔻{T}) where T<:Real = tr(D, H)



"""
    quadraticForm(v::Vector{T}, X::Union(𝕄{T}, ℍ{T}) where T<:RealOrComplex
    quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real

 **alias**: `qf`

 (1) Given a real or complex vector ``v`` and `Hermitian` matrix ``X``,
 compute the quadratic form

 ``v^HXv``,

 where the superscript *H* denotes complex conjugate transpose.

 If ``v`` and ``X`` are real and ``X`` is `Hermitian`, only the
 lower triangular part of ``X`` is used.

 (2) Compute the quadratic form given a real vector ``v``
 and the `LowerTriangular` view ``L`` of a real matrix.

 ## Math

 For ``v`` and ``X`` real and ``X`` symmetric, the quadratic form is

 ``\\sum_i(v_i^2x_{ii})+\\sum_{i>j}(2v_iv_jx_{ij})``.


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
function quadraticForm(v::Vector{T}, X::Union{𝕄{T}, ℍ{T}}) where T<:RealOrComplex
    if T<:Real && X isa(ℍ)
        return quadraticForm(v, 𝕃(X))
    else
        return v'*X*v
    end
end

function quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real
    r=length(v)
    s=0.
    for i=1:r @inbounds s+=(v[i]^2 * L[i, i]) end # Diagonal
    for j=1:r-1, i=j+1:r @inbounds s+=2*v[i]*v[j]*L[i, j] end # Off-diagonal
    return s
end
qufo=quadraticForm


"""
    fidelity(P::ℍ, Q::ℍ)

 Given two positive definte `Hermitian` matrices ``P`` and ``Q``,
 return their *fidelity*:

  ``tr\\big(P^{1/2}QP^{1/2}\\big)^{1/2}.``

  This is used in quantum physics and is related to the
 [Wasserstein](@ref) metric. See for example Bhatia, Jain and Lim (2019b)[🎓](@ref).

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
    (1) fDiagonal(func::Function, X::𝔻, k::Int=0)
    (2) fDiagonal(func::Function, X::𝕃, k::Int=0)
    (3) fDiagonal(func::Function, X::Union{𝕄, ℍ}, k::Int=0)


 **alias**: `𝑓𝔻`

 Applies function `func` element-wise to the elements of the ``k^{th}``
 diagonal of matrix ``X`` (square in all cases but for the 𝕄=`Matrix` argument,
 in which case it may be of dimension *r⋅c*)
 and return a diagonal matrix with these elements.

 See julia [tril(M, k::Integer)](https://bit.ly/2Tbx8o7) function
 for numbering of diagonals.

 Bt default the main diagonal is considered.
 - If the matrix is Diagonal (1) ``k`` must be zero (main diagonal).
 - If the matrix is lower triangular (2) ``k`` cannot be positive.

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
fDiagonal(func::Function, X::𝔻, k::Int=0) = func.(X)

function fDiagonal(func::Function, X::𝕃, k::Int=0)
 if k>0 @error("in function fDiagonal (linearAlgebra.jl): k argument cannot be positive.")
 else return 𝔻(func.(diag(X, k)))
 end
end

fDiagonal(func::Function, X::Union{𝕄, ℍ}, k::Int=0) = 𝔻(func.(diag(X, k)))

𝑓𝔻=fDiagonal



#  -------------------------------
## 5. Unitary functions of matrices
#  -------------------------------
"""
    mgs(T::𝕄, numCol::Int=0)

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
function mgs(T::𝕄, numCol::Int=0)
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

 Given a positive semi-definite `Hermitian` matrix ``S``,
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
    (1) spectralFunctions(P::ℍ, func)
    (2) spectralFunctions(D::𝔻{T}, func) where T<:Real

 (1) This is the *mother function* for all spectral functions of eigenvalues implemented
 in this library, which are:
 - `pow`     (power),
 - `isqrt`   (inverse square root).

 The function `sqr` (square) does not use it, as it can be obtained more
 efficiently by simple multiplication.

 You can use this function if you need another spectral function of eigenvalues
 besides those and those already implemented in the standard package `LinearAlgebra`.
 In general, you won't call it directly.

 `func` is the function that will be applied on the eigenvalues.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).
 It must be a positive definite or positive semi-definite matrix,
 depending on 'func'.

 A special method is provided for real `Diagonal` matrices (2).

!!! note "Nota Bene"
    The function `func` must support the `func.` syntax and therefore
    must be able to apply element-wise to the eigenvalues
    (those include anonymous functions).

 ### Maths

 The definition of spectral functions for a positive definite matrix ``P``
 is at it follows:

 ``f\\big(P\\big)=Uf\\big(Λ\\big)U^H,``

 where ``U`` is the matrix holding in columns the eigenvectors of ``P``,
 ``Λ`` is the matrix holding on diagonal its eigenvalues and ``f`` is
 a function applying element-wise to the eigenvalues.


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
    ispos(F.values, msg="function "*string(func)*": at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return ℍ(F.vectors * 𝔻(func.(F.values)) * F.vectors')
end

spectralFunctions(D::𝔻{T}, func::Function) where T<:Real = func.(D)




"""
    (1) pow(P::ℍ, args...)
    (2) pow(D::𝔻{T}, args...) where T<:Real

 (1) Given a positive semi-definite `Hermitian` matrix ``P``, return the power
 ``P^{r_1}, P^{r_2},...``
 for any number of exponents ``r_1, r_2,...``.
 It returns a tuple of as many elements as arguments passed after ``P``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 ``arg1, arg2,...`` are real numbers.

 A special method is provided for real `Diagonal` matrices (2).

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

pow(D::𝔻{T}, p)  where T<:Real = spectralFunctions(D, x->x^p) # one argument

function pow(P::ℍ, args...)               # several arguments
    (Λ, U) = evd(P)
    ispos(Λ, msg="function Rpow: at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return  (ℍ(U * Λ^p * U') for p in args)
end

function pow(D::𝔻{T}, args...) where T<:Real  # several arguments
    ispos(D, msg="function Rpow: at least one eigenvalue is smaller than the default tolerance")
    return  (D^p for p in args)
end

"""
    (1) invsqrt(P::ℍ)
    (2) invsqrt(D{T}::𝔻) where T<:Real

 Given a positive definite `Hermitian` matrix ``P``,
 compute the inverse of the principal
 square root ``P^{-1/2}``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 A special method is provided for real `Diagonal` matrices (2).

 **See**: [typecasting matrices](@ref).

 **See also**: [`pow`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(ComplexF64, 5);
    Q=invsqrt(P);
    Q*P*Q ≈ I ? println(" ⭐ ") : println(" ⛔ ")

"""
invsqrt(P::ℍ) = spectralFunctions(P, x->1/sqrt(x))

invsqrt(D::𝔻{T}) where T<:Real = spectralFunctions(D, x->1/sqrt(x))



"""
    (1) sqr(P::ℍ)
    (2) sqr(X::Union{𝕄, 𝕃, 𝔻{T}}) where T<:Real

 (1) Given a positive semi-definite `Hermitian` matrix ``P``,
 compute its square ``P^{2}``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 A method is provided also for generic matrices of the `Matrix` type,
 `LowerTriangular` matrices and real `Diagonal` matrices (2). The output
 is of the same type as the input.

 **See also**: [`pow`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(5);
    P²=sqr(P);  # =>  P²=PP
    sqrt(P²)≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
sqr(P::ℍ) = ℍ(P*P)

sqr(X::Union{𝕄, 𝕃, 𝔻{T}}) where T<:Real = X*X


"""
    powerIterations(H::Union{ℍ, 𝕄}, q::Int;
                    <evalues=false, tol::Real=0, maxiter=300, ⍰=false>)

    powerIterations(L::𝕃{T}, q::Int;
                    <evalues=false, tol::Real=0, maxiter=300, ⍰=false)> where T<:Real

 **alias**: `powIter`

 (1) Compute the ``q`` eigenvectors associated to the ``q`` largest (real) eigenvalues
 of real or complex `Hermitian` matrix or `Matrix` ``H`` using the
 [power iterations](https://bit.ly/2JSo0pb) +
 [Gram-Schmidt orthogonalization](https://bit.ly/2YE6zvy) as suggested by Strang.
 The eigenvectors are returned with the same type as the elements of ``H``.

 ``H`` must have real eigenvalues, that is, it must be a symmetric matrix if it is real
 or an Hermitian matrix if it is complex.

 (2) as in (1), but using only the `LowerTriangular` view ``L`` of a matrix.
 This option is available only for real matrices (see below).

 The following are *<optional keyword arguments>*:
 - ``tol`` is the tolerance for the convergence of the power method (see below),
 - ``maxiter`` is the maximum number of iterations allowed for the power method,
 - if ``⍰=true``, the convergence of all iterations will be printed,
 - if ``evalues=true``, return the 4-tuple ``(Λ, U, iterations, covergence)``,
 - if ``evalues=false`` return the 3-tuple ``(U, iterations, covergence)``.


!!! note "Nota Bene"
    Differently from the [`evd`](@ref) function, the eigenvectors and
    eigenvalues are sorted by decreasing order of eigenvalues.

    If ``H`` is real, only its lower triangular part is used like in (2).
    In this case a BLAS routine is used for computing the power iterations.
    See [BLAS routines](@ref).

    ``tol`` defaults to 100 times the square root of `Base.eps` of the type
    of ``H``. This corresponds to requiring the relative convergence criterion
    over two successive iterations to vanish for about half the significant
    digits minus 2.

**See also**: [`mgs`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    # Generate an Hermitian (complex) matrix
    H=randP(ComplexF64, 10);
    # 3 eigenvectors and eigenvalues
    U, iterations, convergence=powIter(H, 3, ⍰=true)
    # all eigenvectors
    Λ, U, iterations, convergence=powIter(H, size(H, 2), evalues=true, ⍰=true);
    U'*U ≈ I && U*Λ*U'≈H ? println(" ⭐ ") : println(" ⛔ ")

    # passing a `Matrix` object
    Λ, U, iterations, convergence=powIter(Matrix(H), 3, evalues=true)

    # passing a `LowerTriangular` object (must be a real matrix in this case)
    L=𝕃(randP(10))
    Λ, U, iterations, convergence=powIter(L, 3, evalues=true)

"""
function powerIterations(H::𝕄, q::Int;
  evalues=false, tol::Real=0, maxiter=300, ⍰=false)

    (n, sqrtn, type) = size(H, 1), √(size(H, 1)), eltype(H)
    tol==0 ? tolerance = √eps(real(type))*1e2 : tolerance = tol
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
        conv = √norm((💡)' * U - I) / sqrtn # relative difference to identity
        (saddlePoint = conv ≈ oldconv)  && @info(msg1, iter)
        (overRun     = iter == maxiter) && @warn(msg2, iter)
        #diverged    = conv > oldconv && @warn(msg3, iter)
        ⍰ && println("iteration: ", iter, "; convergence: ", conv)
        if conv<=tolerance || saddlePoint==true || overRun==true
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
    powerIterations(𝕄(Symmetric(L, :L)), q; evalues=evalues, tol=tol, maxiter=maxiter, ⍰=⍰)

powIter=powerIterations


#  -----------------------------------------------
## 8. Decompositions involving triangular matrices
#  -----------------------------------------------
"""
    (1) choL(P::ℍ)
    (2) choL(D::𝔻{T}) where T<:Real

 (1) Given a real or complex positive definite `Hermitian` matrix ``P``,
 return the *Cholesky lower triangular factor* ``L``
 such that ``LL^H=P``. To obtain ``L^H`` or both ``L`` and ``L^H``, use instead
 julia function [cholesky(P)](https://bit.ly/2u9Hw5P).

 On output, ``L`` is of type [`LowerTriangular`](https://bit.ly/2U511f3).

 (2) For a real `Diagonal` matrix ``D``, return ``D^{1/2}``.

 ## Examples
    using PosDefManifold
    P=randP(5);
    L=choL(P);
    L*L'≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
function choL(P::ℍ)
    choP = cholesky(P)
    return choP.L
end

choL(D::𝔻{T}) where T<:Real = √D
