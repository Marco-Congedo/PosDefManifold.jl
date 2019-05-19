#   Unit linearAlgebra.jl, part of PosDefManifold Package for julia language
#   v 0.0.3 - last update 16th of May 2019
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
#   0. Internal functions
#   1. Utilities
#   2. Matrix normalizations
#   3. Boolean functions of matrices
#   4. Scalar functions of matrices
#   5. Diagonal functions of matrices
#   6. Unitary functions of matrices
#   7. Matrix function of matrices
#   8. Spectral decompositions of positive matrices
#   9. Decompositions involving triangular matrices
# __________________________________________________________________

# -----------------------------------------------------------
# 0. Internal Functions
#    By convention their name begin with underscore char
# -----------------------------------------------------------

# return a vector of ranges partitioning lineraly and
# as much as possible evenly `n` elements in `threads` ranges.
# `threads` is the number of threads to which the ranges are to be
# dispatched. If `threads` is not provided, it is set to the number
# of threads Julia is currently instructed to use.
# For example, for `k`=99
# and `threads`=4, return Array{UnitRange{Int64},1}:[1:25, 26:50, 51:75, 76:99].
# This function is called by threaded function `fVec`
function _partitionLinRange4threads(n::Int, threads::Int=0)
    threads==0 ? thr=nthreads() : thr=threads
    n<thr ? thr = n : nothing
    d = max(round(Int64, n / thr), 1)
    return [(r<thr ? (d*r-d+1:d*r) : (d*thr-d+1:n)) for r=1:thr]
end

# used by function fvec
function _fVec_common(𝐏::AnyMatrixVector;
					  w::Vector=[], ✓w=false, allocs=[])
	threads=nthreads()
	k, n = dim(𝐏, 1), dim(𝐏, 2)
	isempty(w) ? v=[] : v = _getWeights(w, ✓w, k)
	threads==1 || k<threads*4 && @warn fVecMsg threads k
	#allocs==[] ? 𝐐=𝕄Vector([𝕄{type}(undef, n, n) for i=1:threads]) : 𝐐=allocs
	allocs==[] ? 𝐐=𝕄Vector(repeat([𝐏[1]], threads)) : 𝐐=allocs
	return (_partitionLinRange4threads(k, threads), 𝐐, v)
end

#  ------------------------
## 1. Utilities
#  ------------------------
"""
    function typeofMatrix(array::Union{ AnyMatrix,
                                        AnyMatrixVector,
                                        AnyMatrixVector₂ })

 **alias**: `typeofMat`

 Return the type of a matrix, either `Hermitian`,
 `Diagonal`, `LowerTriangular`, or `Matrix`.
 Argument `array` may be a matrix of one of these types, but also one of the following:

 ℍVector, ℍVector₂, 𝔻Vector, 𝔻Vector₂, 𝕃Vector, 𝕃Vector₂, 𝕄Vector, 𝕄Vector₂.

 Those are [Array of Matrices types](@ref).
 See also [aliases](@ref) for the symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.

 Note that this function is different from Julia function
 [typeof](https://docs.julialang.org/en/v1/base/base/#Core.typeof),
 which returns the concrete type (see example below), thus
 cannot be used for [typecasting matrices](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(3) # generate a 3x3 Hermitian matrix
    typeofMatrix(P) # returns `Hermitian`
    typeof(P) # returns `Hermitian{Float64,Array{Float64,2}}`
    # typecast P as a `Matrix` M
    M=Matrix(P)
    # typecast M as a matrix of the same type as P and write the result in A
    A=typeofMatrix(P)(M)

    Pset=randP(3, 4) # generate a set of 4 3x3 Hermitian matrix
    # Pset is an ℍVector type
    typeofMatrix(Pset) # again returns `Hermitian`

"""
typeofMatrix(H::Union{ℍ, ℍVector, ℍVector₂}) = ℍ
typeofMatrix(D::Union{𝔻, 𝔻Vector, 𝔻Vector₂}) = 𝔻
typeofMatrix(L::Union{𝕃, 𝕃Vector, 𝕃Vector₂}) = 𝕃
typeofMatrix(M::Union{𝕄, 𝕄Vector, 𝕄Vector₂}) = 𝕄
typeofMat=typeofMatrix


"""
    function typeofVector(array::Union{ AnyMatrixVector,
                                        AnyMatrixVector₂ })

 **alias**: `typeofVec`

 Return the type of a Vector, either `HermitianVector`,
 `DiagonalVector`, `LowerTriangularVector`, or `MatrixVector`.
 The aliases of those are, respectvely, ℍVector, 𝔻Vector, 𝕃Vector and 𝕄Vector.
 Argument `array` may be a vector of one of these types, but also one of the following:

 ℍVector₂, 𝔻Vector₂, 𝕃Vector₂, 𝕄Vector₂.

 Those are [Array of Matrices types](@ref).
 See also [aliases](@ref) for the symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.

 Note that this function is different from Julia function
 [typeof](https://docs.julialang.org/en/v1/base/base/#Core.typeof)
 only in that it returns the vector type also for
 the ℍVector₂, 𝔻Vector₂, 𝕃Vector₂ and 𝕄Vector₂ types.

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(3, 4) # generate 4 3x3 Hermitian matrix
    typeofMatrix(P) # returns `Array{Hermitian,1}`
    typeof(P) # also returns `Array{Hermitian,1}`

"""
typeofVector(H::Union{ℍVector, ℍVector₂}) = ℍVector
typeofVector(D::Union{𝔻Vector, 𝔻Vector₂}) = 𝔻Vector
typeofVector(L::Union{𝕃Vector, 𝕃Vector₂}) = 𝕃Vector
typeofVector(M::Union{𝕄Vector, 𝕄Vector₂}) = 𝕄Vector
typeofVec=typeofVector


"""
    (1) function dim(X::AnyMatrix, [d])
    (2) function dim(vector::AnyMatrixVector, [d])
    (3) function dim(vector₂::AnyMatrixVector₂, [d])

 (1) ``X`` is a real or complex `Matrix`, `Diagonal`,
 `LowerTriangular` or `Hermitian` matrix.
 Return a 2-tuple containing the dimensions of ``X``,
 which is two times the same dimension for all possible types of ``X``
 with the exception of the `Matrix` type, which can be rectangular.
 Optionally you can specify a dimension (1 or 2)
 to get just the length of that dimension.

 (2) `vector` is an 𝕄Vector, 𝔻Vector, 𝕃Vector or ℍVector type
 (see [AnyMatrixVector type](@ref)).
 Return a 3-tuple containing the number of matrices it holds
 (dimension 1) and their dimensions (dimension 2 and 3).
 Optionally you can specify a dimension (1, 2, or 3)
 to get just the length of that dimension.

 (3) `vector₂` is an 𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂ or ℍVector₂ type
 (see [AnyMatrixVector type](@ref)).
 Return a 4-tuple containing
 - the number of vectors of matrices it holds (dimension 1),
 - a vector holding the number of matrices in each vector of matrices (dimensions 2),
 - the two dimensions of the matrices (dimension 3 and 4).
 Optionally you can specify a dimension (1, 2, 3 or 4)
 to get just the length of that dimension.

 `vector` and `vector₂` are [Array of Matrices types](@ref).
 See also [aliases](@ref) for the symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.

!!! note "Nota Bene"
    If you specify a dimension and this is out of the valid range,
    the function returns zero.

    Both the `vector`(2) and the `vector₂`(3) object are meant to hold
    matrices living in the same manifold, therefore it is assumed
    that all matrices they holds are of the same dimension.
    The dimensions of the matrices are retrived from
    - the first matrix in `vector`(2),
    - the first matrix in the first vector of `vector₂`(3).

 This function replaces Julia [size](https://docs.julialang.org/en/v1/base/arrays/#Base.size)
 function, which cannot be used to retrive dimension for matrix vectors.
 It is not possible to overload the `size` function for matrix vectors
 since this causes problems to other Julia functions.

 ## Examples
    using LinearAlgebra, PosDefManifold
    # (1)
    M=randn(3, 4) # generate a 3x4 `Matrix`
    dim(M) # returns (3, 4)
    dim(M, 1) # returns 3
    dim(M, 2) # returns 4
    dim(M, 3) # out of range: returns 0

    # (2)
    Pset=randP(3, 4) # generate an ℍVector holding 4 3x3 Hermitian matrices
    dim(Pset) # returns (4, 3, 3)
    dim(Pset, 1) # returns 4
    dim(Pset, 2) # returns 3
    dim(Pset, 3) # returns 3

    # (3)
    # Generate a set of 4 random 3x3 SPD matrices
    Pset=randP(3, 4)
    # Generate a set of 40 random 4x4 SPD matrices
    Qset=randP(3, 40)
    A=ℍVector₂([Pset, Qset])
    dim(A) # return (2, [4, 40], 3, 3)
    dim(A, 1) # return 2
    dim(A, 2) # return [4, 40]
    dim(A, 2)[1] # return 4
    dim(A, 3) # return 3
    dim(A, 4) # return 3
    dim(A, 5) # out of range: return 0

    # note: to create an ℍVector₂ object holding k ℍVector objects use
    sets=ℍVector₂(undef, k) # and then fill them
"""
dim(X::AnyMatrix, d::Int) = 1<=d<=2 ? size(X, d) : 0
dim(X::AnyMatrix) = size(X)

function dim(vector::AnyMatrixVector, d::Int) # change to bold X
    if      d==1    return length(vector)
    elseif  2<=d<=3 return size(vector[1], d-1)
    elseif          return 0
    end
end
dim(vector::AnyMatrixVector) = (length(vector), size(vector[1], 1), size(vector[1], 2))

function dim(vector₂::AnyMatrixVector₂, d::Int)
    if      d==1    return length(vector₂)
    elseif  d==2    return collect(length(vec) for vec in vector₂)
    elseif  3<=d<=4 return size(vector₂[1][1], d-2)
    elseif          return 0
    end
end
dim(vector₂::AnyMatrixVector₂) =
    (length(vector₂), collect(length(vec) for vec in vector₂), size(vector₂[1][1], 1), size(vector₂[1][1], 2))

#  ------------------------
## 2. Matrix Normalizations
#  ------------------------

"""
    function det1(X::AnyMatrix; <tol::Real=0>)

 Return the argument matrix ``X`` normalized so as to have *unit determinant*.
 For square positive definite matrices this is the best approximant
 from the set of matrices in the [special linear group](https://bit.ly/2W5jDZ6) -
 see Bhatia and Jain (2014)[🎓].

 ``X`` can be a real or complex `Diagonal`, `LowerTriangular`,
 `Matrix`, or `Hermitian` matrix. (see [AnyMatrix type](@ref))

 If the determinant is not greater to `tol` (which defalts to zero)
 a warning is printed and ``X`` is returned.

!!! note "Nota Bene"
    This function is meant for positive definite matrices.
    Julia may throws an error while computing
    the determinant if the matrix is defective.

 **See** [Julia det function](https://bit.ly/2Y4MnTF).

 **See also**: [`tr1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=det1(P)
    det(Q) # must be 1
    # using a tolerance
    Q=det1(P; tol=1e-12)
"""
function det1(X::AnyMatrix; tol::Real=0)
    tol>=0 ? tolerance=tol : tolerance = 0
    Det = det(X)
    if Det>tolerance X/Det^(1/size(X, 1)) else @warn det1Msg Det tolerance; X end
end
det1Msg="function det1 in LinearAlgebra.jl of PosDefMaifold package: the determinant of the input matrix is not greater than the tolerance."


"""
    tr1(X::AnyMatrix)

 Return the argument matrix ``X`` normalized so as to have *unit trace*.

 ``X`` can be a real or complex `Diagonal`, `LowerTriangular`,
 `Matrix` or `Hermitian` matrix (see [AnyMatrix type](@ref)).

 If the trace is not greater to `tol`
 (which defalts to zero) a warning is printed and ``X`` is returned.

 **See**: [Julia trace function](https://bit.ly/2HoOLiM).

 **See also**: [`tr`](@ref), [`det1`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(5) # generate a random real positive definite matrix 5x5
    Q=tr1(P)
    tr(Q)  # must be 1
    # using a tolerance
    Q=tr1(P; tol=1e-12)

"""
function tr1(X::AnyMatrix; tol::Real=0)
    tol>=0 ? tolerance=tol : tolerance = 0
    trace = tr(X)
    if trace>tolerance X/trace else @warn tr1Msg trace tolerance; X end
end
tr1Msg="function tr1 in LinearAlgebra.jl of PosDefMaifold package: the trace of the input matrix is not greater than the tolerance."


"""
    (1) normalizeCol!(X::𝕄{T}, j::Int)
    (2) normalizeCol!(X::𝕄{T}, j::Int, by::Number)
    (3) normalizeCol!(X::𝕄{T}, range::UnitRange)
    (4) normalizeCol!(X::𝕄{T}, range::UnitRange, by::Number)
                     for all above: where T<:RealOrComplex

 Given a `Matrix` type ``X`` comprised of real or complex elements,
 - (1) normalize the ``j^{th}`` column to unit norm
 - (2) divide the elements of the ``j^{th}`` column by number ``by``
 - (3) normalize the columns in ``range`` to unit norm
 - (4) divide the elements of columns in ``range``  by number ``by``.

 ``by`` is a number of abstract supertype [Number](https://bit.ly/2JwXjGr).
 It should be an integer, real or complex number.
 For efficiency, it should be of the same type as the elements of ``X``.

 ``range`` is a [UnitRange](https://bit.ly/2HSfK5J) type.

 Methods (1) and (3) call the
 [BLAS.nrm2](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.BLAS.nrm2)
 routine for computing the norm of concerned columns.
 See [Threads](@ref).

!!! note "Nota Bene"
    Julia does not allow normalizing the columns of `Hermitian` matrices.
    If you want to call this function for an `Hermitian` matrix see [typecasting matrices](@ref).

 **See** [norm](https://bit.ly/2TaAkR0) and also [randn](https://bit.ly/2I1Vgrg)
 for the example below.

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
function normalizeCol!(X::𝕄{T}, j::Int) where T<:RealOrComplex
    w=colNorm(X, j)
    for i=1:size(X, 1) @inbounds X[i, j]/=w end
end

normalizeCol!(X::𝕄{T}, j::Int, by::Number) where T<:RealOrComplex =
             for i=1:size(X, 1) @inbounds X[i, j]/=by end

function normalizeCol!(X::𝕄{T}, range::UnitRange) where T<:RealOrComplex
    l=range[1]-1
    w=[colNorm(X, j) for j in range]
    for j in range normalizeCol!(X, j, w[j-l]) end
end

normalizeCol!(X::𝕄{T}, range::UnitRange, by::Number) where T<:RealOrComplex =
             for j in range normalizeCol!(X, j, by) end


#  -------------------------------
## 3. Boolean Functions of Matrices
#  -------------------------------

"""
```
    (1) ispos(   λ::Vector{T};
                <tol::Real=0, rev=true, 🔔=true, msg="">) where T<:Real

    (2) ispos(   Λ::𝔻{T};
                <tol::Real=0, rev=true, 🔔=true, msg="">) where T<:Real
```

 Return ``true`` if all numbers in (1) real vector ``λ`` or in (2) real `Diagonal`
 matrix ``Λ`` are not inferior to ``tol``, otherwise return ``false``. This is used,
 for example, in spectral functions to check that all eigenvalues are positive.

!!! note "Nota Bene"
    ``tol`` defaults to the square root of `Base.eps` of the type of ``λ`` (1)
     or ``Λ`` (2). This corresponds to requiring positivity beyond about half of
     the significant digits.

 The following are *<optional keyword arguments>*:

 - If ``rev=true`` the (1) elements in ``λ`` or (2) the diagonal elements in ``Λ`` will be chacked in reverse order.
 This is done for allowing a very fast check when the elements
 are sorted and it is known from where is best to start checking.

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
function ispos(λ::Vector{T}; tol::Real=0, rev=true, 🔔=true, msg="") where T<:Real
    tol==0 ? tolerance = √eps(T) : tolerance = tol
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

ispos(Λ::Diagonal{T}; tol::Real=0, rev=true, 🔔=true, msg="") where T<:Real =
      ispos( diag(Λ); tol=tol, rev=rev, 🔔=🔔, msg=msg)


#  -------------------------------
## 4. Scalar functions of matrices
#  -------------------------------

"""
    (1) colProd(X::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int)
    (2) colProd(X::Union{𝕄{T}, ℍ{T}}, Y::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int)
               for all above: where T<:RealOrComplex

 (1) Given a real or complex `Matrix` or `Hermitian` matrix ``X``,
 return the dot product of the ``j^{th}`` and ``l^{th}`` columns, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*x_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and ``^*`` denotes complex
 conjugate (nothing if the matrix is real).

 (2) Given real or complex `Matrix` or `Hermitian` matrices ``X`` and ``Y``,
 return the dot product of the ``j^{th}`` column of ``X`` and the ``l^{th}`` column
 of ``Y``, defined as,

 ``\\sum_{i=1}^{r} \\big(x_{ij}^*y_{il}\\big), ``

 where ``r`` is the number of rows of ``X`` and of ``Y`` and ``^*`` is as above.

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
colProd(X::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int) where T<:RealOrComplex =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], X[:, l]))

colProd(X::Union{𝕄{T}, ℍ{T}}, Y::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int) where T<:RealOrComplex =
        𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], Y[:, l]))


"""
    colNorm(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex

 Given a real or complex `Matrix` or `Hermitian` matrix ``X``,
 return the Euclidean norm of its ``j^{th}`` column.

 This function calls the
 [BLAS.nrm2](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.BLAS.nrm2)
 routine. See [Threads](@ref).

 **See also**: [`normalizeCol!`](@ref), [`colProd`](@ref), [`sumOfSqr`](@ref).

 ## Examples
    using PosDefManifold
    X=randn(10, 20)
    normOfSecondColumn=colNorm(X, 2)

"""
colNorm(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex =
        BLAS.nrm2(size(X, 1), X[:, j], 1)


"""
    (1) sumOfSqr(A::Array)
    (2) sumOfSqr(H::ℍ{T})
    (3) sumOfSqr(L::𝕃{T})
    (4) sumOfSqr(D::𝔻{T})
    (5) sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, j::Int)
    (6) sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, range::UnitRange)
                 for (1)-(6) above: where T<:RealOrComplex

**alias**: `ss`

 Return
 - (1) the sum of squares of the elements in an array ``A`` of any dimensions.
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

function sumOfSqr(H::ℍ{T}) where T<:RealOrComplex
    r=size(H, 1)
    s=real(T)(0)
    for j=1:size(H, 2)-1
        @inbounds s+=abs2(H[j, j])
        for i=j+1:r @inbounds s+=2*abs2(H[i, j]) end
    end
    @inbounds s+=abs2(H[r, r])
    return s
end

function sumOfSqr(L::𝕃{T}) where T<:RealOrComplex
    s=real(T)(0)
    for j=1:size(L, 2), i=j:size(L, 1)
        @inbounds s+=abs2(L[i, j])
    end
    return s
end

sumOfSqr(D::𝔻{T}) where T<:RealOrComplex = sumOfSqrDiag(D)

sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex = 𝚺(abs2.(X[:, j]))

sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, range::UnitRange) where T<:RealOrComplex =
         𝚺(sumOfSqr(X, j) for j in range)

ss=sumOfSqr


"""
    sumOfSqrDiag(X::AnyMatrix)

 **alias**: `ssd`

 Sum of squares of the diagonal elements in real or complex `Matrix`,
 `Diagonal`, `Hermitian` or `LowerTriangular` matrix ``X``.
 If ``X`` is rectangular (which can be only if it is of the `Matrix` type),
 the main diagonal is considered.

 **See** [AnyMatrix type](@ref)

 **See also**: [`sumOfSqr`](@ref), [`sumOfSqrTril`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    X=randn(10, 20)
    sumDiag2=sumOfSqrDiag(X) # (1)
    sumDiag2=sumOfSqrDiag(𝔻(X)) # (2) 𝔻=LinearAlgebra.Diagonal

"""
sumOfSqrDiag(X::𝕄{T}) where T<:RealOrComplex =
    𝚺(abs2(X[i, i]) for i=1:minimum(size(X)))

sumOfSqrDiag(X::Union{𝔻{T}, ℍ{T}, 𝕃{T}}) where T<:RealOrComplex =
    𝚺(abs2(X[i, i]) for i=1:size(X, 1))

ssd=sumOfSqrDiag

"""
    sumOfSqrTril(X::AnyMatrix, k::Int=0)

**alias**: `sst`

 Given a real or complex `Matrix`, `Diagonal`, `Hermitian` or
 `LowerTriangular` matrix ``X`` (see [AnyMatrix type](@ref)),
 return the sum of squares of the elements
 in its lower triangle up to the ``k^{th}`` underdiagonal.

 `Matrix` ``X`` may be rectangular.

 ``k`` must be in range
 - `1-size(X, 1):c-1` for ``X`` `Matrix`, `Diagonal` or `Hermitian`,
 - `1-size(X, 1):0` for ``X`` `LowerTriangular`.

 For ``X`` `Diagonal` the result is
 - ``0`` if ``k<0``,
 - the sum of the squares of the diagonal elements otherwise.

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
function sumOfSqrTril(X::AnyMatrix, k::Int=0)
    (r, c) = size(X)
    X isa 𝕃 ? range = (1-r:0) : range = (1-r:c-1)
    if k in range
        s=eltype(X)(0)
        for j=1:c, i=max(j-k, 1):r @inbounds s+=abs2(X[i, j]) end
        return s
    else
        @error "in LinearAmgebraInP.sumOfSqrTRil function: argument k is out of bounds"
    end
end
sst=sumOfSqrTril

"""
    (1) tr(P::ℍ{T}, Q::ℍ{T})
    (2) tr(P::ℍ{T}, M::𝕄{T})
    (3) tr(D::𝔻{T}, H::Union{ℍ{T}, 𝕄{T}})
    (4) tr(H::Union{ℍ{T}, 𝕄{T}}, D::𝔻{T})
          for all above: where T<:RealOrComplex

 Given (1) two `Hermitian` positive definite matrix ``P`` and ``Q``,
 return the trace of the product ``PQ``.
 This is real even if ``P`` and ``Q`` are complex.

 ``P`` must always be flagged as `Hermitian`. See [typecasting matrices](@ref).

 In (2) ``Q`` is a `Matrix` object,
 in which case return
 - a real trace if the product ``PQ`` is real or if it has all positive real eigenvalues.
 - a complex trace if the product ``PQ`` is not real and has complex eigenvalues.

 Methods (3) and (4) return the trace of the product ``DH`` or ``HD``,
 where ``D`` is a `Diagonal` matrix and ``H`` an ``Hermitian``
 or ``Matrix`` object. The result is of the same type as the input matrices.

 For all methods all arguments must be of the same type.

 ## Math
 Let ``P`` and ``Q`` be `Hermitian` matrices, using the properties of the trace
 (e.g., the cyclic property and the similarity invariance) you can use this
 function to fast compute the trace of several expressions. For example:

 ``\\textrm{tr}(PQ)=\\textrm{tr}(P^{1/2}QP^{1/2})``

 and

 ``\\textrm{tr}(PQP)=\\textrm{tr}(P^{2}Q)`` (see example below).


 **See**: [trace](https://bit.ly/2HoOLiM).

 **See also**: [`DiagOfProd`](@ref), [`tr1`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
    Q=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
    tr(P, Q) ≈ tr(P*Q) ? println(" ⭐ ") : println(" ⛔ ")
    tr(P, Q) ≈ tr(sqrt(P)*Q*sqrt(P)) ? println(" ⭐ ") : println(" ⛔ ")
    tr(sqr(P), Q) ≈ tr(P*Q*P) ? println(" ⭐ ") : println(" ⛔ ")

"""
tr(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex = real(tr(DiagOfProd(P, Q)))

function tr(P::ℍ{T}, M::𝕄{T}) where T<:RealOrComplex
    λ = [colProd(P, M, i, i) for i=1:size(P, 2)]
    OK=true
    for l in λ
        if imag(l) ≉  0
            OK=false
            break
        end
    end
    if OK return real(𝚺(λ)) else return 𝚺(λ) end
end


function tr(D::𝔻{T}, H::Union{ℍ{T}, 𝕄{T}}) where T<:RealOrComplex
    s=T(0)
    for i=1:size(D, 1) @inbounds s += D[i, i] * H[i, i] end
    return s
end

tr(H::Union{ℍ{T}, 𝕄{T}}, D::𝔻{T}) where T<:RealOrComplex = tr(D, H)



"""
    (1) quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:Real
    (2) quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real
    (3) quadraticForm(v::Vector{T}, X::𝕄{T}, forceLower::Bool=false) where T<:Real
    (4) quadraticForm(v::Vector{S}, X::Union{𝕄{S}, ℍ{S}}) where S<:Complex


 **alias**: `qf`

 (1) Given a real vector ``v`` and a real `Hermitian` matrix ``P``,
 compute the quadratic form

 ``v^TPv``,

 where the superscript *T* denotes transpose,
 using only the lower triangular part of ``P``.

 (2) As in (1), given a real vector ``v``
 and the `LowerTriangular` view ``L`` of a real symmetric matrix.

 (3) As in (1), given a real vector ``v``
 and a real generic `Matrix` ``M``, if `forceLower=true`. If `forceLower=false`,
 the product ``v^TMv`` is evaluated using the whole matrix ``M``.

 (4) Quadratic form ``v^HPv``, where superscript *H* denotes complex conjugate
 and transpose, for a complex vector `v` and  `Matrix` or `Hermitian` matrix.
 The whole matrix is then used.

 ## Math

 For ``v`` and ``X`` real and ``X`` symmetric, the quadratic form is

 ``\\sum_i(v_i^2x_{ii})+\\sum_{i>j}(2v_iv_jx_{ij})``.

 This is used in (1), (2) and (3).


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
quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:Real = qf(v, 𝕃(P))

quadraticForm(v::Vector{T}, X::𝕄{T}, forceLower::Bool=false) where T<: Real = forceLower==true ? qf(v, 𝕃(X)) : v'*X*v

quadraticForm(v::Vector{T}, X::Union{𝕄{T}, ℍ{T}}) where T<:Complex = v'*X*v

function quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real
    r=length(v)
    s=T(0)
    for j=1:r-1
        @inbounds s+=(v[j]^2 * L[j, j])
        for i=j+1:r @inbounds s+=2*v[i]*v[j]*L[i, j]
        end
    end
    @inbounds s+=(v[r]^2 * L[r, r])
    return s
end
qf=quadraticForm


"""
    fidelity(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex

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
function fidelity(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex
    A = √(P)
    return tr(√ℍ(A*Q*A'))
end


#  ---------------------------------
## 5. Diagonal functions of matrices
#  ---------------------------------
"""
    fDiag(func::Function, X::AnyMatrix, k::Int=0)

 **alias**: `𝑓𝔻`

 Applies function `func` element-wise to the elements of the ``k^{th}``
 diagonal of real or complex `Diagonal`, `LowerTriangular`, `Matrix`
 or `Hermitian` matrix ``X`` and return a diagonal matrix with these elements.
 ``X`` must be square in all cases, but for the 𝕄=`Matrix` type argument,
 in which case it may be of dimension *r⋅c*, with *r ≠ c*.

 See julia [tril(M, k::Integer)](https://bit.ly/2Tbx8o7) function
 for numbering of diagonals.

 Bt default the main diagonal is considered.
 - If ``X`` is `Diagonal`, ``k`` must be zero (main diagonal).
 - If ``X`` is `LowerTriangular`, ``k`` cannot be positive.

 Note that if ``X`` is rectangular the dimension of the result depends
 on the size of ``X`` and on the chosen diagonal.
 For example,
 - *r ≠ c* and ``k``=0 (main diagonal), the result will be of dimension min*(r,c)*⋅*min(r,c)*,
 - ``X`` *3⋅4* and ``k=-1``, the result will be *2⋅2*,
 - ``X`` *3⋅4* and ``k=1``, the result will be *3⋅3*, etc.

!!! note "Nota Bene"
    The function `func` must support the `func.` syntax and therefore
    must be able to apply element-wise to the elements of the chosen diagonal
    (this includes [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions-1)).
	If the input matrix is complex, the function `func`
    must be able to support complex arguments.

 **See also**: [`DiagOfProd`](@ref), [`tr`](@ref).

 ## Examples
    using PosDefManifold
    P=randP(5) # use P=randP(ComplexF64, 5) for generating an Hermitian matrix
    D=fDiag(inv, P, -1)   # diagonal matrix with the inverse of the first sub-diagonal of P
    (Λ, U) = evd(P)       # Λ holds the eigenvalues of P, see evd
    Δ=fDiag(log, Λ)       # diagonal matrix with the log of the eigenvalues
    Δ=fDiag(x->x^2, Λ)    # using an anonymous function for the square of the eigenvalues
"""
fDiag(func::Function, X::𝔻{T}, k::Int=0) where T<:RealOrComplex = func.(X)

function fDiag(func::Function, X::𝕃{T}, k::Int=0)  where T<:RealOrComplex
 if k>0 @error("in function fDiag (linearAlgebra.jl): k argument cannot be positive.")
 else return 𝔻(func.(diag(X, k)))
 end
end

fDiag(func::Function, X::Union{𝕄{T}, ℍ{T}}, k::Int=0)  where T<:RealOrComplex =
      𝔻(func.(diag(X, k)))

𝑓𝔻=fDiag


"""
    DiagOfProd(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex

 **alias**: `dop`

 Return the `Diagonal` matrix holding the diagonal of the product ``PQ``
 of two `Hermitian` matrices `P` and `Q`. Only the diagoanl part
 of the product is computed.

 **See also**: [`tr`](@ref), [`fDiag`](@ref).

 ## Examples
    using PosDefManifold, LinearAlgebra
    P, Q=randP(5), randP(5)
    DiagOfProd(P, Q)≈Diagonal(P*Q) ? println("⭐ ") : println("⛔ ")
"""
DiagOfProd(P::ℍ{T}, Q::ℍ{T}) where T<:RealOrComplex =
           𝔻([colProd(P, Q, i, i) for i=1:size(P, 1)])

dop=DiagOfProd


#  -------------------------------
## 6. Unitary functions of matrices
#  -------------------------------
"""
    mgs(X::𝕄{T}, numCol::Int=0) where T<:RealOrComplex

 Modified (stabilized) [Gram-Schmidt orthogonalization](https://bit.ly/2YE6zvy)
 of the columns of square or tall matrix ``X``, which can be comprised of real
 or complex elements.
 The orthogonalized ``X`` is returned by the function. ``X`` is not changed.

 All columns are orthogonalized by default. If instead argument `numCol` is provided,
 then only the first `numCol` columns of ``X`` are orthogonalized.
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
function mgs(X::𝕄{T}, numCol::Int=0) where T<:RealOrComplex
    r=size(X, 1)
    if numCol != 0 && numCol in 2:size(X, 2) c=numCol else c=size(X, 2) end
    U=Matrix{T}(undef, r, c)
    U[:, 1] = X[:, 1]/colNorm(X, 1);
    for i in 2:c
        U[:, i]=X[:, i]
        for j in 1:i-1
            s = colProd(U, j, i) / sumOfSqr(U, j)
            U[:, i] -= s * U[:, j]
        end
        normalizeCol!(U, i)
    end
    return U
end # mgs function


#  ------------------------------
## 7. Matrix function of matrices
#  ------------------------------

"""
	(1) fVec(f::Function, 𝐏::AnyMatrixVector;
		 	<w::Vector=[], ✓w=false, allocs=[]>)

	(2) fVec(f::Function, g::Function, 𝐏::AnyMatrixVector;
			<w::Vector=[], ✓w=false, allocs=[]>)


 Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` matrices
 of the [𝕄Vector type](@ref), [𝔻Vector type](@ref), [𝕃Vector type](@ref) or
 [ℍVector type](@ref) and an optional non-negative real weights vector
 ``w={w_1,...,w_k}``, return expression

 ``(1)\\hspace{6pt}f_{i=1}^{k}(w_iP_i)``,

 or

 ``(2)\\hspace{6pt}f_{i=1}^{k}(w_ig(P_i))``,

 where ``f`` is a linear matrix function iterating over all elements of ``𝐏``,
 typically the `mean` or `sum` function and ``g`` is whatever matrix function
 applying to each matrix ``P_k``, such as `exp`, `log, `sqrt`, etc,
 and [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions-1).

 This function is always **multi-threaded**. It works by partitioning the ``k``
 operations required by the ``f`` function in several groups,
 passing each group to a separate thread and applying the ``f`` function
 again on the intermediate results (that's why ``f`` must be linear).
 This function allows a gain in computational time only when the number of
 matrices (1) and/or their size (2) is high. Use `mean` and `sum` otherwise.
 For the number of threads Julia is instructed to use see [Threads](@ref).

 *<optional keword argument>* `allocs` allows to pass pre-allocated memory
 for holding the intermediate result of each thread.
 Argument `allocs` must be a vector of as many matrices as threads and where
 the matrices have the same dimension as the the matrices in ``𝐏``
 (see the example here below). Using this option is worthwhile only
 if the size of the matrices is very high and/or when `fVec` is to be
 called repeatedly on many vector of matrices, where the matrices
 have always the same size, so that one allocation works for all calls.

 If *<optional keyword argument>* `✓w=true` is passed, the weights are
 normalized so as to sum up to 1, otherwise they are used as they are passed.
 This option is provided to allow calling this function repeatedly without
 normalizing the same weights vector each time. By default `✓w` is false.

!!! note "Nota Bene"
    Contrarily to Julia `mean` and `sum` function (v 1.1.0) the `fVec` function
    returns a matrix of the same type of the matrices in ``𝐏``.

 ## Examples

    using LinearAlgebra, PosDefManifold
    Pset=randP(4, 1000); # generate 1000 positive definite 4x4 matrices
	mean(Pset) # arithmetic mean calling Julia function
	Threads.nthreads() # check that at least two threads are available
	fVec(mean, Pset) multi-threaded arithmetic mean

	inv(mean(inv, Pset)) # Harmonic mean calling Julia function
	inv(fVec(mean, inv, Pset)) # multi-threaded Harmonic mean

	exp(mean(log, Pset)) # log Euclidean mean calling Julia function
	exp(fVec(mean, log, Pset)) # multi-threaded log Euclidean mean

	# notice that Julia `exp` function has changed the type of the result
	# to `Symmetric`. To obtain an `Hermitian` output use
	ℍ(exp(fVec(mean, log, Pset)))

	w=(randn(1000)).^2
	w=w./sum(w)  		# generate normalized random weights

	# weighted arithmetic mean calling Julia function
	sum(Pset[i]*w[i] for i=1:length(w))
	# multi-threaded weighted arithmetic mean
	fVec(sum, Pset, w=w)

	# weighted harmonic mean calling Julia function
	inv(sum(inv(Pset[i])*w[i] for i=1:length(w)))
	# multi-threaded weighted harmonic mean
	inv(fVec(sum, inv, Pset, w=w))

	# pre-allocating memory
	Pset=randP(100, 1000); # generate 1000 positive definite 100x100 matrices
	Qset=MatrixVector(repeat([Pset[1]], nthreads()))
	fVec(mean, log, Pset, allocs=Qset)

	# How much computing time we save ?
	# (example min time obtained with 4 threads & 4 BLAS threads)
	using BenchmarkTools
	# standard Julia function
	@benchmark(mean(log, Pset)) 					# (5.271 s)
	# fVec
	@benchmark(fVec(mean, log, Pset))				# (1.540 s)
"""
function fVec(f::Function, 𝐏::AnyMatrixVector;
			  w::Vector=[], ✓w=false, allocs=[])

	ranges, 𝐐, v = _fVec_common(𝐏; w=w, ✓w=✓w, allocs=allocs)
	l=length(ranges) # number of threads
	if isempty(w)
		@threads for r=1:l 𝐐[r]=f(𝐏[i] for i in ranges[r]) end
	else
		@threads for r=1:l 𝐐[r]=f(v[i]*𝐏[i] for i in ranges[r]) end
	end
    l==1 ? (return 𝐐[1]) : (return typeofMatrix(𝐏)(f(𝐐)))
end

function fVec(f::Function, g::Function, 𝐏::AnyMatrixVector;
			  w::Vector=[], ✓w=false, allocs=[])

	ranges, 𝐐, v = _fVec_common(𝐏; w=w, ✓w=✓w, allocs=allocs)
	l=length(ranges) # number of threads
	if isempty(w)
		@threads for r=1:l 𝐐[r]=f(g(𝐏[i]) for i in ranges[r]) end
	else
		@threads for r=1:l 𝐐[r]=f(v[i]*g(𝐏[i]) for i in ranges[r]) end
	end
    l==1 ? (return 𝐐[1]) : (return typeofMatrix(𝐏)(f(𝐐)))
end
fVecMsg="function fVec of linearAlgebra.jl (PosDefManifold): only one thread is used or the number of matrices (k) is too low for taking advantage of multi-threading"



#  -----------------------------------------------
## 7. Spectral decompositions of positive matrices
#  -----------------------------------------------

"""
    evd(S::ℍ{T}) where T<:RealOrComplex

 Given a positive semi-definite matrix ``S``,
 returns a 2-tuple ``(Λ, U)``, where ``U`` is the matrix holding in columns
 the eigenvectors and ``Λ`` is the matrix holding the eigenvalues on the diagonal.
 This is the output of Julia
 [eigen](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.eigen)
 function in ``UΛU'=S`` form.

 As for the `eigen` function, the eigenvalues and associated
 eigenvectors are sorted by increasing values of eigenvalues.

 ``S`` may be real or complex and must be flagged by Julia as `Hermitian`.
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
function evd(S::ℍ{T}) where T<:RealOrComplex # return tuple (Λ, U)
    F = eigen(S)
    return  𝔻(F.values), F.vectors # 𝔻=LinearAlgebra.Diagonal
end


"""
    (1) spectralFunctions(P::ℍ{T}, func) where T<:RealOrComplex
    (2) spectralFunctions(D::𝔻{S}, func) where S<:Real

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
 depending on `func`.

 A special method is provided for real `Diagonal` matrices (2).

!!! note "Nota Bene"
    The function `func` must support the `func.` syntax and therefore
    must be able to apply element-wise to the eigenvalues
    (those include [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions-1)).

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
function spectralFunctions(P::ℍ{T}, func::Function) where T<:RealOrComplex
    F = eigen(P)
    ispos(F.values, msg="function "*string(func)*": at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return ℍ(F.vectors * 𝔻(func.(F.values)) * F.vectors')
end

spectralFunctions(D::𝔻{T}, func::Function) where T<:Real = func.(D)



"""
    (1) pow(P::ℍ{T}, args...) where T<:RealOrComplex
    (2) pow(D::𝔻{S}, args...) where S<:Real

 (1) Given a positive semi-definite `Hermitian` matrix ``P``, return the power
 ``P^{r_1}, P^{r_2},...``
 for any number of exponents ``r_1, r_2,...``.
 It returns a tuple comprising as many elements as arguments passed after ``P``.

 ``P`` must be flagged as `Hermitian`. See [typecasting matrices](@ref).

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
pow(P::ℍ{T}, p) where T<:RealOrComplex = spectralFunctions(P, x->x^p) # one argument

pow(D::𝔻{T}, p)  where T<:Real = spectralFunctions(D, x->x^p) # one argument

function pow(P::ℍ{T}, args...) where T<:RealOrComplex # several arguments
    (Λ, U) = evd(P)
    ispos(Λ, msg="function pow: at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return  (ℍ(U * Λ^p * U') for p in args)
end

function pow(D::𝔻{T}, args...) where T<:Real  # several arguments
    ispos(D, msg="function pow: at least one eigenvalue is smaller than the default tolerance")
    return  (D^p for p in args)
end

"""
    (1) invsqrt(P{T}::ℍ) where T<:RealOrComplex
    (2) invsqrt(D{S}::𝔻) where S<:Real

 Given a positive definite `Hermitian` matrix ``P``,
 compute the inverse of the principal
 square root ``P^{-1/2}``.

 ``P`` must be flagged as Hermitian. See [typecasting matrices](@ref).

 A special method is provided for real `Diagonal` matrices (2).

 ### Maths

 The principal square root of a positive definite matrix ``P`` is the only
 symmetric (if ``P`` is real) or Hermitian (if ``P`` is complex) square root.
 Its inverse ``P^{-1/2}`` is also named the **whitening** or **sphering**
 matrix since``P^{-1/2}PP^{-1/2}=I``.

 **See**: [typecasting matrices](@ref).

 **See also**: [`pow`](@ref).

 ## Examples
    using LinearAlgebra, PosDefManifold
    P=randP(ComplexF64, 5);
    Q=invsqrt(P);
    Q*P*Q ≈ I ? println(" ⭐ ") : println(" ⛔ ")

"""
invsqrt(P::ℍ{T}) where T<:RealOrComplex = spectralFunctions(P, x->1/sqrt(x))

invsqrt(D::𝔻{T}) where T<:Real = spectralFunctions(D, x->1/sqrt(x))



"""
    (1) sqr(P::ℍ{T}) where T<:RealOrComplex
    (2) sqr(X::Union{𝕄{T}, 𝕃{T}, 𝔻{S}}) where T<:RealOrComplex where S<:Real

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
sqr(P::ℍ{T}) where T<:RealOrComplex = ℍ(P*P)

sqr(X::Union{𝕄{T}, 𝕃{T}, 𝔻{S}}) where T<:RealOrComplex where S<:Real = X*X


"""
    powerIterations(H::Union{ℍ{T}, 𝕄{T}}, q::Int;
    <evalues=false, tol::Real=0, maxiter=300, ⍰=false>)  where T<:RealOrComplex

    powerIterations(L::𝕃{S}, q::Int;
    <evalues=false, tol::Real=0, maxiter=300, ⍰=false)> where S<:Real

 **alias**: `powIter`

 (1) Compute the ``q`` eigenvectors associated to the ``q`` largest (real) eigenvalues
 of real or complex `Hermitian` or `Matrix` ``H`` using the
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

    If ``H`` is `Hermitian` and real, only its lower triangular part is used
    for computing the power iterations, like in (2). In this case the
    [BLAS.symm](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.BLAS.symm)
    routine is used.
    Otherwise the [BLAS.gemm](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.BLAS.gemm)
    routine is used. See [Threads](@ref).

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
function powerIterations(H::𝕄{T}, q::Int;
  evalues=false, tol::Real=0, maxiter=300, ⍰=false) where T<:RealOrComplex

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
        # power iteration U-<H*U of the q vectors of U and their Gram-Schmidt Orth.
        type<:Real ? 💡=mgs(BLAS.symm('L', 'L', H, U)) : 💡=mgs(BLAS.gemm('N', 'N', H, U))
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
        type<:Real ? d=[qf(💡[:, i], H, true) for i=1:q] : d=[qf(💡[:, i], H) for i =1:q]
        return (𝔻(real(d)), 💡, iter, conv)
    end
end


powerIterations(H::ℍ{T}, q::Int;
    evalues=false, tol::Real=0, maxiter=300, ⍰=false) where T<:RealOrComplex =
    powerIterations(Matrix(H), q; evalues=evalues, tol=tol, maxiter=maxiter, ⍰=⍰)

powerIterations(L::𝕃{T}, q::Int;
        evalues=false, tol::Real=0, maxiter=300, ⍰=false) where T<:Real =
    powerIterations(𝕄(L), q; evalues=evalues, tol=tol, maxiter=maxiter, ⍰=⍰)

powIter=powerIterations


#  -----------------------------------------------
## 9. Decompositions involving triangular matrices
#  -----------------------------------------------
"""
    (1) choL(P::ℍ{T}) where T<:RealOrComplex
    (2) choL(D::𝔻{S}) where S<:Real

 (1) Given a real or complex positive definite `Hermitian` matrix ``P``,
 return the *Cholesky lower triangular factor* ``L``
 such that ``LL^H=P``. To obtain ``L^H`` or both ``L`` and ``L^H``, use instead
 julia function [cholesky](https://bit.ly/2u9Hw5P).

 On output, ``L`` is of type [`LowerTriangular`](https://bit.ly/2U511f3).

 (2) For a real `Diagonal` matrix ``D``, return ``D^{1/2}``.

 ## Examples
    using PosDefManifold
    P=randP(5);
    L=choL(P);
    L*L'≈ P ? println(" ⭐ ") : println(" ⛔ ")

"""
function choL(P::ℍ{T}) where T<:RealOrComplex
    choP = cholesky(P)
    return choP.L
end

choL(D::𝔻{T}) where T<:Real = √D
