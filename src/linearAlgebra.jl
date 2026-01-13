#   Unit linearAlgebra.jl, part of PosDefManifold Package for julia language
#
#   MIT License
#   Copyright (c) 2019-25, Marco Congedo, CNRS, Grenobe, France:
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

# return a vector of ranges partitioning linearly and
# as much as possible evenly `n` elements in `threads` ranges.
# `threads` is the number of threads to which the ranges are to be
# dispatched. If `threads` is not provided, it is set to the number
# of threads Julia is currently instructed to use.
# For example, for `k`=99
# and `threads`=4, return Array{UnitRange{Int64},1}:[1:25, 26:50, 51:75, 76:99].
# This function is called by threaded function `fVec`
function _partitionLinRange4threads(n::Int, threads::Int=0)
    threads < 1 ? thr = nthreads() : thr = threads
    #n<thr && throw(ArgumentError("PosDefManifold, internal function `_partitionLinRange4threads`: n must be larger than the number of threads"))
    d = n ÷ thr # integer division
    e = n % thr # reminder
    range = Vector{Any}(undef, thr)
    if e == 0
        for r = 1:thr
            range[r] = d*r-d+1:d*r
        end
    else
        for r = 1:thr-1
            range[r] = d*r-d+1:d*r
        end
        range[thr] = d*(thr-1)+1:n
    end
    return range
end


function _GetThreads(n::Int, callingFunction::String)
    thr = Threads.nthreads()
    n < 1 && throw(ArgumentError("PosDefManifold.jl, internal function `_GetThreads`: `n` must be a positive integer (n=($n))"))
    return n >= thr ? thr : 1
    # return min(n, thr) # this does not work
end

function _GetThreadsAndLinRanges(n::Int, callingFunction::String)
    threads = _GetThreads(n, callingFunction)
    ranges = _partitionLinRange4threads(n, threads)
    return threads, ranges
end

# used by function fvec
function _fVec_common(𝐏::AnyMatrixVector;
    w::Vector=[], ✓w=false, allocs=[])
    threads, ranges = _GetThreadsAndLinRanges(dim(𝐏, 1), "fVec")
    isempty(w) ? v = [] : v = _getWeights(w, ✓w)
    #allocs==[] ? 𝐐=𝕄Vector([𝕄{type}(undef, n, n) for i=1:threads]) : 𝐐=allocs
    #allocs==[] ? 𝐐=𝕄Vector(repeat([zeros(eltype(𝐏[1]), size(𝐏[1]))], threads)) : 𝐐=allocs
    allocs == [] ? 𝐐 = 𝕄Vector(repeat([similar(𝐏[1])], threads)) : 𝐐 = allocs
    return (threads, ranges, 𝐐, v)
end

#  ------------------------
## 1. Utilities
#  ------------------------
"""
    function typeofMatrix(
	array::Union{AnyMatrix, AnyMatrixVector, AnyMatrixVector₂})

**alias**: `typeofMat`

Return the type of a matrix, either `Hermitian`,
`Diagonal`, `LowerTriangular`, or `Matrix`.
Argument `array` may be a matrix of one of these types, but also one of the following:

`ℍVector`, `ℍVector₂`, `𝔻Vector`, `𝔻Vector₂`, `𝕃Vector`, `𝕃Vector₂`,
`𝕄Vector`, `𝕄Vector₂`.

Those are [Array of Matrices types](@ref).
See also [aliases](@ref) for the symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.

Note that this function is different from Julia function
[typeof](https://docs.julialang.org/en/v1/base/base/#Core.typeof),
which returns the concrete type (see example below), thus
cannot be used for [typecasting matrices](@ref).

**Examples**

```julia
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
```

"""
typeofMatrix(H::Union{ℍ, ℍVector, ℍVector₂}) = ℍ
typeofMatrix(D::Union{𝔻, 𝔻Vector, 𝔻Vector₂}) = 𝔻
typeofMatrix(L::Union{𝕃, 𝕃Vector, 𝕃Vector₂}) = 𝕃
typeofMatrix(M::Union{𝕄, 𝕄Vector, 𝕄Vector₂}) = 𝕄
typeofMat=typeofMatrix


"""
    function typeofVector(
	array::Union{AnyMatrix, AnyMatrixVector, AnyMatrixVector₂})

 **alias**: `typeofVec`

Return the type of a Vector, either `HermitianVector`,
`DiagonalVector`, `LowerTriangularVector`, or `MatrixVector`.
The aliases of those are, respectvely, `ℍVector`, `𝔻Vector`, `𝕃Vector` and
`𝕄Vector`.
Argument `array` may be a vector of one of these types, but also one of the
following:

`ℍ`, `𝔻`, `𝕃` and `𝕄`, `ℍVector₂`, `𝔻Vector₂`, `𝕃Vector₂`, `𝕄Vector₂`.

See [aliases](@ref) for the symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.
The last four are [Array of Matrices types](@ref).

Note that this function is different from Julia function
[typeof](https://docs.julialang.org/en/v1/base/base/#Core.typeof)
only in that it returns the vector type also if `array`
is not of the `ℍVector`, `𝔻Vector`, `𝕃Vector` or
`𝕄Vector` type.

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(3, 4) # generate 4 3x3 Hermitian matrix
typeofMatrix(P) # returns `Array{Hermitian,1}`
typeof(P) # also returns `Array{Hermitian,1}`

typeofMatrix(P[1]) # returns `Array{Hermitian,1}`
typeof(P[1]) # returns `Hermitian{Float64,Array{Float64,2}}`
```

"""
typeofVector(H::Union{ℍ, ℍVector, ℍVector₂}) = ℍVector
typeofVector(D::Union{𝔻, 𝔻Vector, 𝔻Vector₂}) = 𝔻Vector
typeofVector(L::Union{𝕃, 𝕃Vector, 𝕃Vector₂}) = 𝕃Vector
typeofVector(M::Union{𝕄, 𝕄Vector, 𝕄Vector₂}) = 𝕄Vector
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

 **Examples**
```julia
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
```

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


"""
```
function remove(X::Union{Vector, Matrix}, what::Union{Int, Vector{Int}};
				dims=1)
```

Remove one or more elements from a vector or one or more
columns or rows from a matrix.

If `X` is a Matrix, `dims`=1 (default) remove rows,
`dims`=2 remove columns.

If `X` is a Vector, `dims` has no effect.

The second argument is either an integer or a vector of integers.

 **Examples**
```julia
a=randn(5)
b=remove(a, 2) # remove second element
b=remove(a, collect(1:3)) # remove rows 1 to 3

A=randn(3, 3)
B=remove(A, 2) # remove second row
B=remove(A, 2; dims=2) # remove second column

A=randn(5, 5)
B=remove(A, collect(1:2:5)) # remove rows 1, 3 and 5
C=remove(A, [1, 4]) # remove rows 1 and 4

# remove columns 2, 3, 8, 9, 10
A=randn(10, 10)
B=remove(A, [collect(2:3); collect(8:10)]; dims=2)

# remove every other sample (decimation by a factor of 2)
A=randn(10, 10)
B=remove(A, collect(1:2:size(A, 1)); dims=1)
```

"""
function remove(X::Union{Vector, Matrix}, what::Union{Int, Vector{Int}}; dims=1)
    1<dims<2 && throw(ArgumentError("function `remove`: the `dims` keyword argument must be 1 or 2"))
    di = X isa Vector ? 1 : dims
    d = size(X, di)
    mi, ma = minimum(what), maximum(what)
    (1≤mi≤d && 1≤ma≤d) || throw(ArgumentError("function `remove`: the second argument must holds elements comprised in between 1 and $d. Check also the `dims` keyword"))
    b = filter(what isa Int ? x->x≠what : x->x∉what, 1:d)
    return X isa Vector ? X[b] : X[di==1 ? b : 1:end, di==2 ? b : 1:end]
end


"""
    function isSquare(X::Matrix)=size(X, 1)==size(X, 2)

Return true if matrix `X` is square, false otherwise.
"""
isSquare(X::Matrix)=size(X, 1)==size(X, 2)


#  ------------------------
## 2. Matrix Normalizations
#  ------------------------

"""
    function det1(X::AnyMatrix; <tol::Real=0.>)

Return the argument matrix ``X`` normalized so as to have *unit determinant*.
For square positive definite matrices this is the best approximant
from the set of matrices in the [special linear group](https://bit.ly/2W5jDZ6) -
see Bhatia and Jain (2014)[🎓](@ref).

``X`` can be a real or complex `Diagonal`, `LowerTriangular`,
`Matrix`, or `Hermitian` matrix. (see [AnyMatrix type](@ref))

If the determinant is not greater than `tol` (which defalts to zero)
a warning is printed and ``X`` is returned.

!!! note "Nota Bene"
    This function is meant for positive definite matrices.
    Julia may throws an error while computing
    the determinant if the matrix is defective.

**See** [Julia det function](https://bit.ly/2Y4MnTF).

**See also**: [`tr1`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(5) # generate a random real positive definite matrix 5x5
Q=det1(P)
det(Q) # must be 1
# using a tolerance
Q=det1(P; tol=1e-12)
```

"""
function det1(X::AnyMatrix; tol::Real=0.)
    tol>=0. ? tolerance=tol : tolerance = 0.
    Det = det(X)
    if Det>tolerance X/Det^(1/size(X, 1)) else @warn det1Msg Det tolerance; X end
end
det1Msg="function det1 in LinearAlgebra.jl of PosDefMaifold package: the determinant of the input matrix is not greater than the tolerance."


"""
    tr1(X::AnyMatrix; tol::Real=0.)

Return the argument matrix ``X`` normalized so as to have *unit trace*.

``X`` can be a real or complex `Diagonal`, `LowerTriangular`,
`Matrix` or `Hermitian` matrix (see [AnyMatrix type](@ref)).
Its trace must be real. If the absolute value of its imaginary part
is greater than `tol` (which defalts to zero) a warning is printed
and ``X`` is returned.
Also, if the trace is not greater than `tol`
a warning is printed and ``X`` is returned.

**See**: [Julia trace function](https://bit.ly/2HoOLiM).

**See also**: [`tr`](@ref), [`det1`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold

P=randP(5) # generate a random real positive definite matrix 5x5
Q=tr1(P)
tr(Q)  # must be 1
# using a tolerance
Q=tr1(P; tol=1e-12)

Pc=randP(ComplexF64, 5) # generate a random real positive definite matrix 5x5
Qc=tr1(Pc)
tr(Qc)  # must be 1
```

"""
function tr1(X::AnyMatrix; tol::Real=0.)
    tol >= 0. ? tolerance = tol : tolerance = 0.
    trace = tr(X)
    imagtr = imag(trace)
    if abs(imagtr) > tolerance
        @warn tr1Msg2 imagtr tolerance
        return X
    else
        trace = real(trace)
    end
    if trace > tolerance
        return X / trace
    else
        @warn tr1Msg1 trace tolerance
        return X
    end
end
tr1Msg1="function tr1 in LinearAlgebra.jl of PosDefMaifold package: the trace of the input matrix is not greater than the tolerance."
tr1Msg2="function tr1 in LinearAlgebra.jl of PosDefMaifold package: the imaginary part of the trace of the input matrix is greater than the tolerance."


"""
    nearestPosDef(X::Union{𝔻, 𝕄}; tol::Real=0.)

Return the nearest symmetric/Hermitian positive semi-definite matrix
of a diagonal or of an arbitary square matrix `X` according to the
Frobenius norm.
If the eigenvalues of the symmetric part of `X` are all non-negative,
the result is positive definite and will be flagged as `Hermitian`,
otherwise it is positive semi-definite and will not be flagged.
The nearest matrix is given by

``(Y+H)/2``

where

``Y=(X+X^H)/2``

is the symmetric part of ``X``, and ``H`` is the symmetric polar factor
of ``Y``. See Higham(1988)[🎓](@ref) for details and for the way it is computed.

**See also**: [`det1`](@ref), [`procrustes`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
X=randn(5, 5) # generate an arbitrary 5x5 matrix
S=nearestPosDef(X)

P=randP(5) # generate a random real positive definite 5x5 matrix
S=nearestPosDef(Matrix(P)) # typecasting an Hermitian matrix as a `Matrix`
# Since P is a positive definite matrix S must be equal to P
S ≈ P ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function nearestPosDef(D::𝔻; tol::Real=0.)
    tol >= 0. ? tolerance = tol : tolerance = 0.
    return 𝔻([D[i, i] >= tolerance ? D[i, i] : 0 for i = 1:size(D, 1)])
end

function nearestPosDef(X::𝕄; tol::Real=0.)
    size(X, 1) == size(X, 2) || throw(ArgumentError("PosDefManifold.jl, function nearestPosDef: the input matrix must be square"))
    tol >= 0. ? tolerance = tol : tolerance = 0.
    F = eigen((X + X') / 2)
    λispos = ispos(F.values; 🔔=false, rev=false)
    D = λispos ? 𝔻(F.values) : nearestPosDef(𝔻(F.values), tol=tolerance)
    return λispos ? ℍ(F.vectors * D * F.vectors') : (F.vectors * D * F.vectors')
end


"""
    nearestOrthogonal(X::AnyMatrix)

**alias**: `nearestOrth`

Return the nearest orthogonal matrix
of a square `Hermitian`, `LowerTriangular`, `Diagonal` or generic `Matrix` `X`
(see [AnyMatrix type](@ref)).
This is given by

``UV^H``,

where

``\\textrm(SVD)=UΛV^H``.

If `X` is `Diagonal`, return `X`.

**See also**: [`nearestPosDef`](@ref), [`procrustes`](@ref).

**Examples**
```julia
using PosDefManifold
U=nearestOrth(randn(5, 5))
```

"""
function nearestOrthogonal(X::AnyMatrix)
    size(X, 1) == size(X, 2) || throw(ArgumentError("PosDefManifold.jl, function nearestOrthogonal: the input matrix must be square"))
    if X isa Diagonal
        return X
    else
        sv = svd(X)
        return BLAS.gemm('N', 'N', sv.U, sv.Vt) # sv.U * sv.Vt
    end
end
nearestOrth = nearestOrthogonal



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

**Examples**
```julia
using PosDefManifold
X=randn(10, 20)
normalizeCol!(X, 2)                  # (1) normalize columns 2
normalizeCol!(X, 2, 10.0)            # (2) divide columns 2 by 10.0
normalizeCol!(X, 2:4)                # (3) normalize columns 2 to 4
X=randn(ComplexF64, 10, 20)
normalizeCol!(X, 3)                  # (1) normalize columns 3
normalizeCol!(X, 3:6, (2.0 + 0.5im)) # (4) divide columns 3 to 5 by (2.0 + 0.5im)
```

"""
function normalizeCol!(X::𝕄{T}, j::Int) where T<:RealOrComplex
    w = colNorm(X, j)
    for i = 1:size(X, 1)
        @inbounds X[i, j] /= w
    end
end

normalizeCol!(X::𝕄{T}, j::Int, by::Number) where T<:RealOrComplex =
    for i = 1:size(X, 1)
        @inbounds X[i, j] /= by
    end

function normalizeCol!(X::𝕄{T}, range::UnitRange) where T<:RealOrComplex
    l = range[1] - 1
    w = [colNorm(X, j) for j in range]
    for j in range
        normalizeCol!(X, j, w[j-l])
    end
end

normalizeCol!(X::𝕄{T}, range::UnitRange, by::Number) where T<:RealOrComplex =
    for j in range
        normalizeCol!(X, j, by)
    end


#  -------------------------------
## 3. Boolean Functions of Matrices
#  -------------------------------

"""
```
    (1) ispos(λ::Vector{T};
	<
	tol::Real=0,
	rev=true,
	🔔=true,
	msg="">)

    (2) ispos(Λ::𝔻{T};
	< same optional keyword arguments as in (1) > )

	for all above: where T<:Real
```

Return ``true`` if all numbers in (1) real vector ``λ`` or in (2) real `Diagonal`
matrix ``Λ`` are not inferior to ``tol``, otherwise return ``false``. This is used,
for example, in spectral functions to check that all eigenvalues are positive.

!!! note "Nota Bene"
    ``tol`` defaults to the square root of `Base.eps` of the type of ``λ`` (1)
     or ``Λ`` (2). This corresponds to requiring positivity beyond about half of
     the significant digits.

The following are *<optional keyword arguments>*:

- If ``rev=true`` the (1) elements in ``λ`` or (2) the diagonal elements in ``Λ`` will be chacked in reverse order. This is done for allowing a very fast check when the elements are sorted and it is known from where is best to start checking.

If the result is ``false``:
- if ``🔔=true`` a bell character will be printed. In most systems this will ring a bell on the computer.
- if string ``msg`` is provided, a warning will print ``msg`` followed by:
"at position *pos*", where *pos* is the position where the
first non-positive element has been found.

**Examples**
```julia
using PosDefManifold
a=[1, 0, 2, 8]
ispos(a, msg="non-positive element found")

# it will print:
# ┌ Warning: non-positive element found at position 2
# └ @ [here julie will point to the line of code issuing the warning]
```
"""
function ispos(λ::Vector{T};
        tol::Real=0,
        rev=true,
        🔔=true,
        msg="") where T<:Real

    tol == 0 ? tolerance = √eps(T) : tolerance = tol
    rev ? iterations = (length(λ):-1:1) : iterations = (1:length(λ))
    for i in iterations
        if λ[i] < tolerance
            🔔 && print('\a') # print('\a') sounds a bell
            length(msg) > 0 && @warn("function ispos(linearAlgebra.jl) " * msg * " at position $i")
            return false
            break
        end
    end
    return true
end

ispos(Λ::Diagonal{T};
        tol::Real=0,
        rev=true,
        🔔=true,
        msg="") where T<:Real =
    ispos(diag(Λ); tol=tol, rev=rev, 🔔=🔔, msg=msg)


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

**Examples**
```julia
using PosDefManifold
X=randn(10, 20)
p=colProd(X, 1, 3)
Y=randn(10, 30)
q=colProd(X, Y, 2, 25)
```

"""
function colProd(X::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int) where T<:RealOrComplex
    s = zero(T)
    @simd for i=1:size(X, 1)
        @inbounds s += conj(X[i, j])*X[i, l]
    end
    s
end  # old version: 𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], X[:, l]))

function colProd(X::Union{𝕄{T}, ℍ{T}}, Y::Union{𝕄{T}, ℍ{T}}, j::Int, l::Int) where T<:RealOrComplex
    s = zero(T)
    @simd for i=1:min(size(X, 1), size(Y, 1))
        @inbounds s += conj(X[i, j])*Y[i, l]
    end
    s
end #old version: 𝚺(conj(x1)*x2 for (x1, x2) in zip(X[:, j], Y[:, l]))


"""
    colNorm(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex

Given a real or complex `Matrix` or `Hermitian` matrix ``X``,
return the Euclidean norm of its ``j^{th}`` column.

This function calls the
[BLAS.nrm2](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.BLAS.nrm2)
routine. See [Threads](@ref).

**See also**: [`normalizeCol!`](@ref), [`colProd`](@ref), [`sumOfSqr`](@ref).

**Examples**
```julia
using PosDefManifold
X=randn(10, 20)
normOfSecondColumn=colNorm(X, 2)
```

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

**Examples**
```julia
using PosDefManifold
X=randn(10, 20)
sum2=sumOfSqr(X)        # (1) sum of squares of all elements
sum2=sumOfSqr(X, 1)     # (2) sum of squares of elements in column 1
sum2=sumOfSqr(X, 2:4)   # (3) sum of squares of elements in column 2 to 4
```
"""
sumOfSqr(A::Array) = 𝚺(abs2(a) for a in A)

function sumOfSqr(H::ℍ{T}) where T<:RealOrComplex
    r=size(H, 1)
    s=real(T)(0)
    for j=1:size(H, 2)-1
        @inbounds s+=abs2(H[j, j])
        @simd for i=j+1:r 
            @inbounds s+=2*abs2(H[i, j]) 
        end
    end
    @inbounds s+=abs2(H[r, r])
    return s
end

function sumOfSqr(L::𝕃{T}) where T<:RealOrComplex
    s=real(T)(0)
    for j=1:size(L, 2) 
        @simd for i=j:size(L, 1)
            @inbounds s+=abs2(L[i, j])
        end
    end
    return s
end

sumOfSqr(D::𝔻{T}) where T<:RealOrComplex = sumOfSqrDiag(D)

sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, j::Int) where T<:RealOrComplex = 𝚺(abs2.(X[:, j]))

sumOfSqr(X::Union{𝕄{T}, ℍ{T}}, range::UnitRange) where T<:RealOrComplex =
         𝚺(sumOfSqr(X, j) for j in range)

ss=sumOfSqr


"""
    (1) sumOfSqrDiff(H1, H2::ℍ{T})
    (2) sumOfSqrDiff(L1, L2::𝕃{T})
    (3) sumOfSqrDiff(D1, D2::𝔻{T})
    for (1)-(3) above: where T<:RealOrComplex

**alias**: `ssdiff`

Return the sum of the squares of the elements of the difference of the arguments, i.e., of
- (1) ``H_1-H_2`` (Hermitian matrices) 
- (2) ``L_1-L_2`` (lower triangular matrices) 
- (3) ``D_1-D_2`` (Diagonal matrices) 

All methods support real and complex matrices.

**Examples**
```julia
using PosDefManifold
H1, H2 = randn(10, 20), randn(10, 20)
sum2=sumOfSqrDiff(H1, H2) 
```
"""
function sumOfSqrDiff(H1, H2::ℍ{T}) where T<:RealOrComplex
    r=size(H1, 1)
    s=real(T)(0)
    for j=1:size(H1, 2)-1
        @inbounds s+=abs2(H1[j, j]-H2[j, j])
        for i=j+1:r @inbounds s+=2*abs2(H1[i, j]-H2[i, j]) end
    end
    @inbounds s+=abs2(H1[r, r]-H2[r, r])
    return s
end

function sumOfSqrDiff(L1, L2::𝕃{T}) where T<:RealOrComplex
    r=size(L1, 1)
    s=real(T)(0)
    for j=1:size(L1, 2)-1
        @inbounds s+=abs2(L1[j, j]-L2[j, j])
        for i=j+1:r @inbounds s+=abs2(L1[i, j]-L2[i, j]) end
    end
    @inbounds s+=abs2(L1[r, r]-L2[r, r])
    return s
end

function sumOfSqrDiff(D1, D2::𝔻{T}) where T<:RealOrComplex
    s=real(T)(0)
    @simd for j=1:size(D1, 2)
        @inbounds s+=abs2(D1[j, j]-D2[j, j])
    end
    return s
end

ssdiff=sumOfSqrDiff


"""
    sumOfSqrDiag(X::AnyMatrix)

**alias**: `ssd`

Sum of squares of the diagonal elements in real or complex `Matrix`,
`Diagonal`, `Hermitian` or `LowerTriangular` matrix ``X``.
If ``X`` is rectangular (which can be only if it is of the `Matrix` type),
the main diagonal is considered.

**See** [AnyMatrix type](@ref)

**See also**: [`sumOfSqr`](@ref), [`sumOfSqrTril`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
X=randn(10, 20)
sumDiag2=sumOfSqrDiag(X) # (1)
sumDiag2=sumOfSqrDiag(𝔻(X)) # (2)
# 𝔻=LinearAlgebra.Diagonal is declated in the main module
```
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

**Examples**
```julia
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
```
"""
function sumOfSqrTril(X::AnyMatrix, k::Int=0)
    (r, c) = size(X)
    X isa 𝕃 ? range = (1-r:0) : range = (1-r:c-1)
    if k in range
        s=eltype(X)(0)
        for j=1:c
            @simd for i=max(j-k, 1):r 
                @inbounds s+=abs2(X[i, j]) 
            end
        end
        return real(s)
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


**See**: [tr](https://bit.ly/2HoOLiM).

**See also**: [`DiagOfProd`](@ref), [`tr1`](@ref).

**Examples**
```julia
using PosDefManifold
P=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
Q=randP(ComplexF64, 5) # generate a random complex positive definite matrix 5x5
tr(P, Q) ≈ tr(P*Q) ? println(" ⭐ ") : println(" ⛔ ")
tr(P, Q) ≈ tr(sqrt(P)*Q*sqrt(P)) ? println(" ⭐ ") : println(" ⛔ ")
tr(sqr(P), Q) ≈ tr(P*Q*P) ? println(" ⭐ ") : println(" ⛔ ")
```

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
    if OK 
        return real(𝚺(λ)) 
    else 
        return 𝚺(λ) 
    end
end


function tr(D::𝔻{T}, H::Union{ℍ{T}, 𝕄{T}}) where T<:RealOrComplex
    s=T(0)
    @simd for i=1:size(D, 1) 
        @inbounds s += D[i, i] * H[i, i] 
    end
    return s
end

tr(H::Union{ℍ{T}, 𝕄{T}}, D::𝔻{T}) where T<:RealOrComplex = tr(D, H)



"""
    (1) quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:Real
    (2) quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real
    (3) quadraticForm(v::Vector{T}, X::𝕄{T}, forceLower::Bool=false) where T<:Real
    (4) quadraticForm(v::Vector{S}, X::Union{𝕄{S}, ℍ{S}, 𝕃{S}}) where S<:Complex

**alias**: `qf`

(1) Given a real vector ``v`` and a real `Hermitian` matrix ``P``,
compute the quadratic form

``v^TPv``,

where the superscript *T* denotes transpose.
It uses only the lower triangular part of ``P``.

(2) As in (1), given a real vector ``v``
and a `LowerTriangular` matrix ``L``.

(3) As in (1), given a real vector ``v``
and a real generic `Matrix` ``M``, if `forceLower=true`. If `forceLower=false`,
the product ``v^TMv`` is evaluated instead using the whole matrix ``M``.

(4) Quadratic form ``v^HPv``, where superscript *H* denotes complex conjugate
and transpose, for a complex vector `v` and a complex `Matrix`, `LowerTrianglar`
or `Hermitian` matrix.
The whole matrix is used.

## Math

For ``v`` and ``X`` real and ``X`` symmetric, the quadratic form is

``\\sum_i(v_i^2x_{ii})+\\sum_{i>j}(2v_iv_jx_{ij})``.

For ``L`` lower triangular is

``\\sum_i(v_i^2x_{ii})+\\sum_{i>j}(v_iv_jx_{ij})``.

These formula are used in methods (1), (2) and (3).


**Examples**
```julia
using PosDefManifold
P=randP(5) # generate a random real positive definite matrix 5x5
v=randn(5)
q1=quadraticForm(v, P) # or q1=qf(v, P)
q2=v'*P*v
q1 ≈ q2 ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function quadraticForm(v::Vector{T}, P::ℍ{T}) where T<:Real
    r = length(v)
    s = T(0)
    for j = 1:r-1
        @inbounds s += (v[j]^2 * P[j, j])
        @simd for i = j+1:r
            @inbounds s += 2 * v[i] * v[j] * P[i, j]
        end
    end
    @inbounds s += (v[r]^2 * P[r, r])
    return s
end

function quadraticForm(v::Vector{T}, X::𝕄{T}, forceLower::Bool=false) where T<:Real
    if forceLower
        r = length(v)
        s = T(0)
        for j = 1:r-1
            @inbounds s += (v[j]^2 * X[j, j])
            @simd for i = j+1:r
                @inbounds s += 2 * v[i] * v[j] * X[i, j]
            end
        end
        @inbounds s += (v[r]^2 * X[r, r])
        return s
    else
        return v' * X * v
    end
end

quadraticForm(v::Vector{T}, X::Union{𝕄{T}, ℍ{T}, 𝕃{T}}) where T<:Complex = v'*X*v

function quadraticForm(v::Vector{T}, L::𝕃{T}) where T<:Real
    r=length(v)
    s=T(0)
    for j=1:r-1
        @inbounds s+=(v[j]^2 * L[j, j])
        @simd for i=j+1:r @inbounds s+=2*v[i]*v[j]*L[i, j] end
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

**Examples**
```julia
using PosDefManifold
P=randP(5);
Q=randP(5);
f=fidelity(P, Q)
```

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
- If ``X`` is `Diagonal`, ``k`` is set automatically to zero (main diagonal).
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

**Examples**
```julia
using PosDefManifold
P=randP(5) # use P=randP(ComplexF64, 5) for generating an Hermitian matrix

# diagonal matrix with the inverse of the first sub-diagonal of P
D=fDiag(inv, P, -1)

(Λ, U) = evd(P) # Λ holds the eigenvalues of P, see evd

# diagonal matrix with the log of the eigenvalues
Δ=fDiag(log, Λ)

# using an anonymous function for the square of the eigenvalues
Δ=fDiag(x->x^2, Λ)
```

"""
fDiag(func::Function, X::𝔻{T}, k::Int=0) where T<:RealOrComplex = 𝔻(func.(diag(X)))

function fDiag(func::Function, X::𝕃{T}, k::Int=0) where T<:RealOrComplex
    if k > 0
        @error("in function fDiag (linearAlgebra.jl): k argument cannot be positive.")
    else
        return 𝔻(func.(diag(X, k)))
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

**Examples**
```julia
using PosDefManifold, LinearAlgebra
P, Q=randP(5), randP(5)
DiagOfProd(P, Q)≈Diagonal(P*Q) ? println("⭐ ") : println("⛔ ")
```

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

**Examples**
```julia
using LinearAlgebra, PosDefManifold
X=randn(10, 10);
U=mgs(X)        # result is 10⋅10
U=mgs(X, 3)     # result is 10⋅3
U'*U ≈ I ? println(" ⭐ ") : println(" ⛔ ")
# julia undertands also:
U'U ≈ I ? println(" ⭐ ") : println(" ⛔ ")
```

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
```
	(1) fVec(f::Function, 𝐏::AnyMatrixVector;
	<
	w::Vector=[],
	✓w=false,
	allocs=[])
	>

	(2) fVec(f::Function, g::Function, 𝐏::AnyMatrixVector;
	< same optional keyword arguments in (1) >)
```

Given a 1d array ``𝐏={P_1,...,P_k}`` of ``k`` matrices
of the [𝕄Vector type](@ref), [𝔻Vector type](@ref), [𝕃Vector type](@ref) or
[ℍVector type](@ref) and an optional non-negative real weights vector
``w={w_1,...,w_k}``, return expression

``(1)\\hspace{6pt}f_{i=1}^{k}(w_iP_i)``,

or

``(2)\\hspace{6pt}f_{i=1}^{k}(w_ig(P_i))``,

where ``f`` is either the `mean` or the `sum` standard julia functions
and ``g`` is whatever matrix function
applying to each matrix ``P_k``, such as `exp`, `log`, `sqrt`, etc,
and [anonymous functions](https://docs.julialang.org/en/v1/manual/functions/#man-anonymous-functions-1).

This function is **multi-threaded**. It works by partitioning the ``k``
operations required by the ``f`` function in several groups,
passing each group to a separate thread and combining the result
of the intermediate operations.
This function allows a gain in computational time only when the number of
matrices (1) and/or their size (2) is high. Use `mean` and `sum` otherwise.
The maximal gain is obtained when the number of matrices in `𝐏` is an exact
multiple of the number of threads Julia is instructed to use.
For this latter, see [Threads](@ref).

!!! note "Nota Bene"
    Contrarily to Julia `mean` and `sum` function (v 1.1.0) the `fVec` function
    returns a matrix of the same type of the matrices in ``𝐏``.

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

**See also** [`typeofMatrix`](@ref), [`typeofVector`](@ref), [`dim`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
Pset=randP(4, 1000); # generate 1000 positive definite 4x4 matrices
mean(Pset) # arithmetic mean calling Julia function
Threads.nthreads() # check how many threads are available
fVec(mean, Pset) # multi-threaded arithmetic mean

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
Qset=MatrixVector(repeat([similar(Pset[1])], Threads.nthreads()))
fVec(mean, log, Pset, allocs=Qset)

# How much computing time we save ?
# (example min time obtained with 4 threads & 4 BLAS threads)
using BenchmarkTools
# standard Julia function
@benchmark(mean(log, Pset)) 					# (5.271 s)
# fVec
@benchmark(fVec(mean, log, Pset))				# (1.540 s)
```

"""
function fVec(f::Function, g::Function, 𝐏::AnyMatrixVector;
        w::Vector=[],
        ✓w=false,
        allocs=[])
    f≠mean && f≠sum && begin
	    @error "the `f` argument of the fVec function must be `mean` or `sum`"
		return
	end
	threads, ranges, 𝐐, v = _fVec_common(𝐏; w=w, ✓w=✓w, allocs=allocs)
	#println("treads ", threads)
	#println("ranges ", ranges)

	if isempty(w)
		g==identity ? (@threads for r=1:threads 𝐐[r]=sum(𝐏[i] for i in ranges[r]) end) :
 			(@threads for r=1:threads 𝐐[r]=sum(g(𝐏[i]) for i in ranges[r]) end)
	else
		g==identity ? (@threads for r=1:threads 𝐐[r]=sum(v[i]*𝐏[i] for i in ranges[r]) end) :
			(@threads for r=1:threads 𝐐[r]=sum(v[i]*g(𝐏[i]) for i in ranges[r]) end)
	end
	if f==sum
    	return threads==1 ? typeofMatrix(𝐏)(𝐐[1]) : typeofMatrix(𝐏)(sum(𝐐))
    else
		invn=inv(length(𝐏))
		return threads==1 ? typeofMatrix(𝐏)(𝐐[1]*invn) : typeofMatrix(𝐏)(sum(𝐐)*invn)
	end
end


fVec(f::Function, 𝐏::AnyMatrixVector;
        w::Vector=[],
        ✓w=false,
        allocs=[]) =
    fVec(f, identity, 𝐏; w=w, ✓w=✓w, allocs=allocs)

@doc raw"""
	(1) congruence(B::AnyMatrix, P::AnyMatrix, matrixType)
	(2) congruence(B::AnyMatrix, 𝐏::AnyMatrixVector, matrixVectorType)
	(3) congruence(B::AnyMatrix, 𝑷::AnyMatrixVector₂, matrixVector₂Type)
	(4) congruence(𝐁::AnyMatrixVector, 𝑷::AnyMatrixVector₂, matrixVector₂Type)

**alias**: `cong`

(1) Return the congruent transformation

``BPB^H``,

for ``B`` and ``P`` any combination of `Hermitian`, `LowerTriangular`,
`Diagonal` or general `Matrix` type.

The result is of the `matrixType` argument, which must be provided and
must be one of these four abstract type (not an instance of them).
See [aliases](@ref) for shortening these type using symbols `ℍ`, `𝔻`, `𝕃` and `𝕄`.

(2) Return a vector of matrices holding the congruent transformations

``BP_kB^H``,

for all ``k`` matrices in ``𝐏={P_1,...,P_k}``, for ``B`` and ``𝐏``
any combination of matrix type `Hermitian`, `LowerTriangular`,
`Diagonal` or `Matrix` (``B``) and vector of matrices type `ℍVector`, `𝔻Vector`,
`𝕃Vector` and `𝕄Vector` (``𝐏``). See [Array of Matrices types](@ref).

The result is a vector of matrices of the `matrixVectorType` argument,
which must be provided and must be one of the following abstract types:
`ℍVector`, `𝔻Vector`, `𝕃Vector` or `𝕄Vector`
(and not an instance of these types).

(3) Return a vector of vector of matrices holding the
congruent transformations

``BP_{mk}B^H``,

for all ``m`` vectors of ``k[m]`` vectors of matrices in ``𝑷``,
for ``𝐁`` and ``𝑷`` any combination of matrix type `Hermitian`,
`LowerTriangular`, `Diagonal` or `Matrix` (``𝐁``) and vector of matrices type
`ℍVector₂`, `𝔻Vector₂`,
`𝕃Vector₂` and `𝕄Vector₂` (``𝑷``). See [Array of Matrices types](@ref).

The result is a vector of vector of matrices of the `matrixVector₂Type`
argument, which must be provided and must be one of the following
abstract types: `ℍVector₂`, `𝔻Vector₂`, `𝕃Vector₂` or `𝕄Vector₂`
(and not an instance of these types).

(4) Return a vector of vector of matrices holding the
congruent transformations

``B_iP_{ij}B_j^H``, for ``i,j∈[1,...,m]``.

for ``𝑩`` holding ``m`` matrices and ``𝑷`` holding ``m`` vectors
holding ``m`` matrices each.
Note that, differently from method (3), here the vectors of ``𝑷``
are all of the same length and this is eaxctly the length of ``𝑩``.
``𝑩`` and ``𝑷`` may be any combination of matrix vector type `ℍVector`,
`𝔻Vector`, `𝕃Vector` and `𝕄Vector` (``𝑩``) and vector of matrices type
`ℍVector₂`, `𝔻Vector₂`, `𝕃Vector₂` and `𝕄Vector₂` (``𝑷``).
See [Array of Matrices types](@ref).

Note that this function computes the following algebraic expression:

``\begin{pmatrix} B_1 & \hspace{0.01cm} & 0 \\ \hspace{0.01cm} & \ddots & \hspace{0.01cm} \\ 0 & \hspace{0.01cm} & B_m \end{pmatrix}
\begin{pmatrix} P_{11} & \cdots & P_{1m} \\ \vdots & \ddots & \vdots \\ P_{m1} & \cdots & P_{mm} \end{pmatrix}
\begin{pmatrix}B_1^T & \hspace{0.01cm} & 0 \\ \hspace{0.01cm} & \ddots & \hspace{0.01cm} \\ 0 & \hspace{0.01cm} & B_m^T\end{pmatrix}``


The result is a vector of vector of matrices of the `matrixVector₂Type`
argument, which must be provided and must be one of the following
abstract types: `ℍVector₂`, `𝔻Vector₂`, `𝕃Vector₂` or `𝕄Vector₂`
(and not an instance of these types).

When you pass it to this function, make sure to typecast ``𝑩``
as an `ℍVector`, `𝔻Vector`, `𝕃Vector` or `𝕄Vector` type if it is not
already created as one of these types. See the example here below
and [typecasting matrices](@ref).

Method (2), (3) and (4) are **multi-threaded**. See [Threads](@ref).

!!! note "Nota Bene"
    Types `ℍ`, `𝔻`, `𝕃` or `𝕄` are actually constructors, thus they may
    modify the result of the congruence(s). This greatly expand the
    possibilities of this function, but it is your responsibility to
    pick the right argument `matrixType` in (1), `matrixVectorType` in (2) and
	`matrixVector₂Type` in (3)-(4).
    For example, in (1) if ``B`` and ``P`` are `Hermitian`,
    calling `cong(B, P, 𝔻)` will actually
    return the diagonal part of ``B*P*B'`` and calling `cong(B, P, 𝕃)` will
    actually return its lower triangular part. The full congruence can
    be obtained as an `Hermitian` matrix by `cong(B, P, ℍ)` and as a generic
    matrix object by `cong(B, P, 𝕄)`.

**See also** [`typeofMatrix`](@ref), [`typeofVector`](@ref), [`dim`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold

# (1)
P=randP(3) # generate a 3x3 positive matrix
M=randn(3, 3)
C=cong(M, P, ℍ) # equivalent to C=ℍ(M*P*M')

# (2)
Pset=randP(4, 100); # generate 100 positive definite 4x4 matrices
M=randn(4, 4)
Qset=cong(M, Pset, ℍVector) # = [M*Pset_1*M',...,M*Pset_k*M'] as an ℍVector type

# recenter the matrices in Pset to their Fisher mean:
Qset=cong(invsqrt(mean(Fisher, Pset)), Pset, ℍVector)

# as a check, the Fisher mean of Qset is now the identity
mean(Fisher, Qset)≈I ? println("⭐") : println("⛔")

# (3)
Pset1=randP(4, 10); # generate 10 positive definite 4x4 matrices
Pset2=randP(4, 8);
Pset=ℍVector₂([Pset1, Pset2]);
M=randn(4, 4)
Qset=cong(M, Pset, MatrixVector₂)
Qset[1][1]≈M*Pset[1][1]*M' ? println("⭐") : println("⛔")
Qset[1][5]≈M*Pset[1][5]*M' ? println("⭐") : println("⛔")
Qset[2][1]≈M*Pset[2][1]*M' ? println("⭐") : println("⛔")
Qset[2][4]≈M*Pset[2][4]*M' ? println("⭐") : println("⛔")

# (4)
Pset1=randP(4, 2); # generate 2 positive definite 4x4 matrices
Pset2=randP(4, 2);
Pset=ℍVector₂([Pset1, Pset2]);
U=𝕄Vector([randU(4), randU(4)])
Qset=cong(U, Pset, MatrixVector₂)
Qset[1][1]≈U[1]*Pset[1][1]*U[1]' ? println("⭐") : println("⛔")
Qset[1][2]≈U[1]*Pset[1][2]*U[2]' ? println("⭐") : println("⛔")
Qset[2][1]≈U[2]*Pset[2][1]*U[1]' ? println("⭐") : println("⛔")
Qset[2][2]≈U[2]*Pset[2][2]*U[2]' ? println("⭐") : println("⛔")
```

"""
congruence(B::AnyMatrix, P::AnyMatrix, matrixType) = matrixType(B*P*B')

function congruence(B::AnyMatrix, 𝐏::AnyMatrixVector, matrixVectorType)
    k, 𝕋 = dim(𝐏, 1), typeofMat(matrixVectorType(undef, 0))

    threads = _GetThreads(k, "congruence")
    if threads == 1
        return matrixVectorType([congruence(B, P, 𝕋) for P in 𝐏])
    else
        𝐐 = matrixVectorType(undef, k)
        @threads for i = 1:k
            𝐐[i] = congruence(B, 𝐏[i], 𝕋)
        end
        return 𝐐
    end
end


function congruence(B::AnyMatrix, 𝑷::AnyMatrixVector₂, matrixVector₂Type)
    m, k, 𝕋 = dim(𝑷, 1), dim(𝑷, 2), typeofVec(matrixVector₂Type(undef, 0)) #NB: k is a vector
    threads = _GetThreads(m, "congruence")
    𝓠 = matrixVector₂Type(undef, m)

    if threads == 1
        for i = 1:m
            𝓠[i] = congruence(B, 𝑷[i], 𝕋)
        end
    else
        @threads for i = 1:m
            𝓠[i] = congruence(B, 𝑷[i], 𝕋)
        end
    end
    return 𝓠
end


function congruence(𝐁::AnyMatrixVector, 𝑷::AnyMatrixVector₂, matrixVector₂Type)
    m, k, dummy = dim(𝑷, 1), dim(𝑷, 2), matrixVector₂Type(undef, 0) #NB: k is a vector
    𝕊, 𝕋 = typeofMat(dummy), typeofVec(dummy)
    threads = _GetThreads(sum(m * k), "congruence")
    𝓠 = matrixVector₂Type(undef, m)
    for i = 1:m
        𝓠[i] = 𝕋(undef, k[i])
    end

    if threads == 1
        for i = 1:m, j = 1:k[i]
            𝓠[i][j] = 𝕊(𝐁[i] * 𝑷[i][j] * 𝐁[j]')
        end
    else
        @threads for i = 1:m
            @threads for j = 1:k[i]
                𝓠[i][j] = 𝕊(𝐁[i] * 𝑷[i][j] * 𝐁[j]')
            end
        end
    end
    return 𝓠
end


cong=congruence



#  -----------------------------------------------
## 7. Spectral decompositions of positive matrices
#  -----------------------------------------------

"""
    evd(S::Union{𝕄{T}, ℍ{T}}) where T<:RealOrComplex

Given a positive semi-definite matrix ``S``,
returns a 2-tuple ``(Λ, U)``, where ``U`` is the matrix holding in columns
the eigenvectors and ``Λ`` is the matrix holding the eigenvalues on the diagonal.
This is the output of Julia
[eigen](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.eigen)
function in ``UΛU'=S`` form.

As for the `eigen` function, the eigenvalues and associated
eigenvectors are sorted by increasing values of eigenvalues.

``S`` may be real or complex and may be flagged by Julia as `Hermitian`
(in this case **PosDefManifold** assumes it is positive definite).

See [typecasting matrices](@ref).

**See also**: [`spectralFunctions`](@ref).

**Examples**
```julia
using PosDefManifold
A=randn(3, 3);
S=A+A';
Λ, U=evd(S); # which is equivalent to (Λ, U)=evd(P)
(U*Λ*U') ≈ S ? println(" ⭐ ") : println(" ⛔ ")
# => UΛU'=S, UΛ=SU, ΛU'=U'S
```

"""
function evd(S::Union{𝕄{T}, ℍ{T}}) where T<:RealOrComplex # return tuple (Λ, U)
    F = eigen(S)
    return  𝔻(F.values), F.vectors # 𝔻=LinearAlgebra.Diagonal
end



"""
```
frf(P::ℍ{T}) where T<:RealOrComplex
```
Full-rank factorization of `Hermitian` matrix `P`.
It is given by

``F=UD^{1/2}``,

where

``\\textrm{EVD}(P)=UDU^{H}``

is the eigenvalue-eigenvector decomposition of `P`. It verifies

``FF^H=P``,

thus ``F^{-1}`` is a whitening matrix.

**See also**: [`invfrf`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(3)
F = frf(P)
F*F'≈P ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function frf(P::ℍ{T}) where T<:RealOrComplex
   #size(P, 1)≠size(A, 2) && throw(ArgumentError(📌*", frf function: input matrix must be square"))
   λ, U=eigen(P)
   return U*Diagonal(sqrt.(λ))
end


"""
```
invfrf(P::ℍ{T}) where T<:RealOrComplex
```
Inverse of the full-rank factorization of `Hermitian` matrix `P`.
It is given by

``F=D^{-1/2}U^H``,

where

``\\textrm{EVD}(P)=UDU^{H}``

is the eigenvalue-eigenvector decomposition of `P`. It verifies

``FPF^H=I``,

thus ``F`` is a whitening matrix.

**See also**: [`frf`](@ref).

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(3)
F = invfrf(P)
F*P*F'≈I ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function invfrf(P::ℍ{T}) where T<:RealOrComplex
   #size(P, 1)≠size(A, 2) && throw(ArgumentError(📌*", frf function: input matrix must be square"))
   Λ, U=evd(P)
   return invsqrt(Λ)*U'
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

**Examples**
```julia
using LinearAlgebra, PosDefManifold
n=5
P=randP(n) # P=randP(ComplexF64, 5) to generate an Hermitian complex matrix
noise=0.1;
Q=spectralFunctions(P, x->x+noise) # add white noise to the eigenvalues
tr(Q)-tr(P) ≈ noise*n ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function spectralFunctions(P::ℍ{T}, func::Function) where T<:RealOrComplex
    F = eigen(P)
    ispos(F.values, msg="function "*string(func)*": at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return ℍ((F.vectors * 𝔻(func.(F.values))) * F.vectors')
end

spectralFunctions(D::𝔻{T}, func::Function) where T<:Real = fDiag(func, D)



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

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(5);     # use P=randP(ComplexF64, 5) for generating an Hermitian matrix
Q=pow(P, 0.5);            # =>  QQ=P
Q, W=pow(P, 0.5, -0.5);
W*P*W ≈ I ? println(" ⭐ ") : println(" ⛔ ")
Q*Q ≈ P ? println(" ⭐ ") : println(" ⛔ ")
R, S=pow(P, 0.3, 0.7);
R*S ≈ P ? println(" ⭐ ") : println(" ⛔ ")
```

"""
pow(P::ℍ{T}, p) where T<:RealOrComplex = spectralFunctions(P, x->x^p) # one argument

pow(D::𝔻{T}, p)  where T<:Real = D^p # one argument

function pow(P::ℍ{T}, args...) where T<:RealOrComplex # several arguments
    (Λ, U) = evd(P)
    ispos(Λ, msg="function pow: at least one eigenvalue is smaller than the default tolerance")
    # optimize by computing only the upper trinagular part
    return  (ℍ((U * Λ^p) * U') for p in args)
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

**Examples**
```julia
using LinearAlgebra, PosDefManifold
P=randP(ComplexF64, 5);
Q=invsqrt(P);
Q*P*Q ≈ I ? println(" ⭐ ") : println(" ⛔ ")
```

"""
invsqrt(P::ℍ{T}) where T<:RealOrComplex = spectralFunctions(P, x->1/√x)

invsqrt(D::𝔻{T}) where T<:Real = spectralFunctions(D, x->1/√x)



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

**Examples**
```julia
using PosDefManifold
P=randP(5);
P²=sqr(P);  # =>  P²=PP
sqrt(P²)≈ P ? println(" ⭐ ") : println(" ⛔ ")
```

"""
sqr(P::ℍ{T}) where T<:RealOrComplex = ℍ(P*P)

sqr(X::Union{𝕄{T}, 𝕃{T}, 𝔻{S}}) where T<:RealOrComplex where S<:Real = X*X


"""
    powerIterations(H::Union{ℍ{T}, 𝕄{T}}, q::Int;
    <
	evalues=false,
	tol::Real=0,
	maxiter::Int=300,
	verbose=false>) where T<:RealOrComplex

    powerIterations(L::𝕃{S}, q::Int;
    < same optional keyword arguments in (1)>) where S<:Real

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
- `tol` is the tolerance for the convergence of the power method (see below),
- `maxiter` is the maximum number of iterations allowed for the power method,
- if `verbose`=true, the convergence of all iterations will be printed,
- if `evalues`=true, return the 4-tuple ``(Λ, U, iterations, covergence)``,
- if `evalues`=false return the 3-tuple ``(U, iterations, covergence)``.

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

**Examples**
```julia
using LinearAlgebra, PosDefManifold
# Generate an Hermitian (complex) matrix
H=randP(ComplexF64, 10);
# 3 eigenvectors and eigenvalues
U, iterations, convergence=powIter(H, 3, verbose=true)
# all eigenvectors
Λ, U, iterations, convergence=powIter(H, size(H, 2), evalues=true, verbose=true);
U'*U ≈ I && U*Λ*U'≈H ? println(" ⭐ ") : println(" ⛔ ")

# passing a `Matrix` object
Λ, U, iterations, convergence=powIter(Matrix(H), 3, evalues=true)

# passing a `LowerTriangular` object (must be a real matrix in this case)
L=𝕃(randP(10))
Λ, U, iterations, convergence=powIter(L, 3, evalues=true)
```

"""
function powerIterations(H::𝕄{T}, q::Int;
        evalues=false,
        tol::Real=0,
        maxiter::Int=300,
        verbose=false) where T<:RealOrComplex

    (n, sqrtn, type) = size(H, 1), √(size(H, 1)), eltype(H)
    tol==0 ? tolerance = √eps(real(type))*1e2 : tolerance = tol
    msg1 = "Power Iterations reached a saddle point at:"
    msg2 = "Power Iterations reached the max number of iterations at:"
    U = randn(type, n, q) # initialization
    normalizeCol!(U, 1:q)
    💡 = similar(U) # 💡 is the poweriteration matrix
    (iter, conv, oldconv) = 1, 0., maxpos
    verbose && @info("Running Power Iterations...")
    while true
        # power iteration U-<H*U of the q vectors of U and their Gram-Schmidt Orth.
        type<:Real ? 💡=mgs(BLAS.symm('L', 'L', H, U)) : 💡=mgs(BLAS.gemm('N', 'N', H, U))
        conv = √norm((💡)' * U - I) / sqrtn # relative difference to identity
        (saddlePoint = conv ≈ oldconv)  && @info(msg1, iter)
        (overRun     = iter == maxiter) && @warn(msg2, iter)
        #diverged    = conv > oldconv && @warn(msg3, iter)
        verbose && println("iteration: ", iter, "; convergence: ", conv)
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
        evalues=false,
        tol::Real=0,
        maxiter::Int=300,
        verbose=false) where T<:RealOrComplex =
    powIter(Matrix(H), q; evalues=evalues, tol=tol, maxiter=maxiter, verbose=verbose)

powerIterations(L::𝕃{T}, q::Int;
        evalues=false,
        tol::Real=0,
        maxiter::Int=300,
        verbose=false) where T<:Real =
    powIter(𝕄(L), q; evalues=evalues, tol=tol, maxiter=maxiter, verbose=verbose)

powIter = powerIterations


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

**See also**: [`choInv`](@ref).

**Examples**
```julia
using PosDefManifold
P=randP(5);
L=choL(P);
L*L'≈ P ? println(" ⭐ ") : println(" ⛔ ")
```

"""
function choL(P::ℍ{T}) where T<:RealOrComplex
    choP = cholesky(P)
    return choP.L
end

choL(D::𝔻{T}) where T<:Real = √D


"""
    choInv(P::AbstractArray{T};
		kind::Symbol = :LLt, tol::Real = √eps(T)) where T<:RealOrComplex

For a real or complex positive definite matrix ``P``, let ``P=LL^H`` be its
*Cholesky decomposition* and ``P=L_1DL_1^H`` the related *LDLt* decomposition.
In the above, ``L`` is a lower triangular matrix, ``D`` a positive-definite
diagonal matrix and ``L_1`` a unit lower triangular matrix.
Return:
- if `kind`is `:LLt` (default), the 2-tuple ``L``, ``L^{-H}``
- if `kind`is `:LDLt`, the 3-tuple ``L_1``, ``D``, ``L_1^{-H}``.

Those are obtained in one pass and for small matrices this is faster
then calling Julia's
[chelosky](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.cholesky)
function and inverting the lower factor unless you set

    BLAS.set_num_threads(1).

Input matrix `P` may be of type `Matrix` or `Hermitian`. Since only the
lower triangle is used, `P` may also be a `LowerTriangular` view of a
positive definite matrix.
If `P` is real, it can also be of the `Symmetric` type.

The algorithm is a *multiplicative Gaussian elimination*.
If run completely, in input matrix `P` there will be the Identity at the end.

**Notes:**
Output ``L^{-H}`` is an inverse square root (whitening matrix) of ``P``,
since ``L^{-1}PL^{-H}=I``. It therefore yields the inversion of ``P`` as
``P^{-1}=L^{-H}L^{-1}``. It is the fastest whitening matrix to be computed,
however it yields poor numerical precision, especially for large matrices.

The following relations holds:
- ``L=PL^{-H}``
- ``L^{H}=L^{-1}P``
- ``L^{-H}=P^{-1}L``
- ``L^{-1}=L^{H}P^{-1}``.

We also have
- ``L^{H}L=L^{-1}P^{2}L^{-H}=UPU^H``, with ``U`` orthogonal (see below) and
- ``L^{-1}L^{-H}=L^{H}P^{-2}L=UP^{-1}U^H``.

``LL^{H}`` and ``L^{H}L`` are unitarily similar, that is,

``ULL^{H}U^H=L^{H}L``,

where ``U=L^{-1}P^{1/2}``, with ``P^{1/2}=H`` the *principal*
(unique symmetric) square root of ``P``. This is seen writing
``PP^{-1}=HHL^{-H}L^{-1}``; multiplying both sides on the left by ``L^{-1}``
and on the right by ``L`` we obtain

``L^{-1}PP^{-1}L=L^{-1}HHL^{-H}=I=(L^{-1}H)(L^{-1}H)^H``

and since ``L^{-1}H`` is square it must be unitary.

From these expressions we have
- ``H=LU=U^HL^H``
- ``L=HU^H``
- ``H^{-1}=U^HL^{-1}``
- ``L^{-1}=UH^{-1}``.

``U`` is the *polar factor* of ``L^{H}``, *i.e.*, ``L^{H}=UH``,
since ``LL^{H}=HU^HUH^H=H^2=P``.

From ``L^{H}L=UCU^H`` we have ``L^{H}LU=UC=ULL^{H}`` and from
``U=L^{-1}H`` we have ``L=HU^H``.

**See also**: [`choInv!`](@ref), [`choL`](@ref).

**Examples**
```julia
using PosDefManifold
n, t = 800, 6000
etol = 1e-9
Z=randn(t, n)
Y=Z'*Z
Yi=inv(Y)

A, B=choInv!(copy(Y))
norm(A*A'-Y)/√n < etol ? println(" ⭐ ") : println(" ⛔ ")
norm(B*B'-Yi)/√n < etol ? println(" ⭐ ") : println(" ⛔ ")

A, D, B=choInv!(copy(Y); kind=:LDLt)
norm(Y-A*D*A')/√n < etol ? println(" ⭐ ") : println(" ⛔ ")
norm(Yi-B*inv(D)*B')/√n < etol ? println(" ⭐ ") : println(" ⛔ ")

# repeat the test for complex matrices
Z=randn(ComplexF64, t, n)
Y=Z'*Z
Yi=inv(Y)

A, B=choInv!(copy(Y))
norm(A*A'-Y)/√n < etol ? println(" ⭐ ") : println(" ⛔ ")
norm(B*B'-Yi)/√n < etol ? println(" ⭐ ") : println(" ⛔ ")

A, D, B=choInv!(copy(Y); kind=:LDLt)
norm(Y-A*D*A')/√n < etol ? println(" ⭐ ") : println(" ⛔ ")
norm(Yi-B*inv(D)*B')/√n < etol ? println(" ⭐ ") : println(" ⛔ ")
```

"""
choInv(P::AbstractArray{T};
	   kind::Symbol = :LLt, tol::Real = √eps(T)) where T<:Real =
    choInv!(P isa Hermitian || P isa Symmetric ? copy(Matrix(P)) : copy(P);
			kind=kind, tol=tol)

choInv(P::AbstractArray{T};
	   kind::Symbol=:LLt, tol::Real = √eps(real(T))) where T<:Complex =
	choInv!(P isa Hermitian ? copy(Matrix(P)) : copy(P);
			kind=kind, tol=tol)

"""
    choInv!(P::AbstractArray{T};
		kind::Symbol = :LLt, tol::Real = √eps(T)) where T<:RealOrComplex

The same thing as [`choInv`](@ref), but destroys the input matrix.
This function does not require copying the input matrix,
thus it is slightly faster.
"""
function choInv!(P::AbstractArray{T};
			  	 kind::Symbol = :LLt, tol::Real = √eps(T)) where T<:Real

    P isa Matrix || P isa LowerTriangular || throw(ArgumentError("function choInv!: input matrix must be of the Matrix or LowerTriangular type. Call `choinv` instead"))
    n = size(P, 1)
    L₁ = kind == :LDLt ? UnitLowerTriangular(zeros(T, n, n)) : LowerTriangular(Matrix{T}(I, n, n))
    U₁⁻¹ = kind == :LDLt ? UnitUpperTriangular(zeros(T, n, n)) : UpperTriangular(Matrix{T}(I, n, n))

    
    for j = 1:n-1
        @inbounds P[j, j] < tol && throw(LinearAlgebra.PosDefException(1))
        for i = j+1:n
            @inbounds θ = P[i, j] / -P[j, j]
            @simd for k = i:n
                @inbounds P[k, i] += θ * P[k, j]
            end # update A and write D
            @inbounds L₁[i, j] = -θ
            @simd for k = 1:j-1
                @inbounds U₁⁻¹[k, i] += θ * U₁⁻¹[k, j]
            end
            @inbounds U₁⁻¹[j, i] = θ
        end
    end

    kind == :LDLt ? (return L₁, Diagonal(P), U₁⁻¹) : begin
        D = sqrt.(Diagonal(P))
        return L₁ * D, U₁⁻¹ * inv(D)
    end
end


function choInv!(P::AbstractArray{T};
    kind::Symbol=:LLt, tol::Real=√eps(real(T))) where T<:Complex
    P isa Matrix || P isa LowerTriangular || throw(ArgumentError("function choInv!: input matrix must be of the Matrix or LowerTriangular type Call `choInv` instead"))
    n = size(P, 1)
    L₁ = kind == :LDLt ? UnitLowerTriangular(zeros(T, n, n)) : LowerTriangular(Matrix{T}(I, n, n))
    U₁⁻¹ = kind == :LDLt ? UnitUpperTriangular(zeros(T, n, n)) : UpperTriangular(Matrix{T}(I, n, n))

    for j = 1:n-1
        @inbounds abs2(P[j, j]) < tol && throw(LinearAlgebra.PosDefException(1))
        for i = j+1:n
            @inbounds θ = conj(P[i, j] / -P[j, j])
            @simd for k = i:n
                @inbounds P[k, i] += θ * P[k, j]
            end # update A and write D
            @inbounds L₁[i, j] = conj(-θ)
            @simd for k = 1:j-1
                @inbounds U₁⁻¹[k, i] += θ * U₁⁻¹[k, j]
            end
            @inbounds U₁⁻¹[j, i] = θ
        end
    end

    kind == :LDLt ? (return L₁, Diagonal(P), U₁⁻¹) : begin
        D = sqrt.(Diagonal(P))
        return L₁ * D, U₁⁻¹ * inv(D)
    end
end
