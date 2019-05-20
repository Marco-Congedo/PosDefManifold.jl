# linearAlgebra.jl

 This unit contains linear algebra functions useful in relation to the Riemannian
 geometry of the manifold of *Symmetric Positive Definite (SPD)* or
 *Hermitian Positive Definite (HPD)* matrices. In Julia those are `Hermitian` matrices, see [typecasting matrices](@ref).

 In general they take a matrix as input (some may take other arrays as input) and are divided in eight categories depending on what kind of functions thay are and what they give as output:

| Category  | Output |
|:---------- |:----------- |
| 1. [Utilities](@ref)  | - - - |
| 2. [Matrix normalizations](@ref) | matrix |
| 3. [Boolean functions of matrices](@ref) | matrix |
| 4. [Scalar functions of matrices](@ref) | scalar |
| 5. [Diagonal functions of matrices](@ref) | diagonal matrix |
| 6. [Unitary functions of matrices](@ref) | orthogonal/unitary matrix |
| 7. [Matrix function of matrices](@ref) | matrix |
| 8. [Spectral decompositions of positive matrices](@ref) | spectral function of input|
| 9. [Decompositions involving triangular matrices](@ref) | triangular matrix |

⋅

## Utilities

| Function   | Description |
|:---------- |:----------- |
| [`typeofMatrix`](@ref), `typeofMat` | Return the type of the matrix argument |
| [`typeofVector`](@ref), `typeofVec` | Return the type of the matrix vector argument |
| [`dim`](@ref)| length of the dimansions of matrices and vectors of matrices |


⋅

```@docs
typeofMatrix
typeofVector
dim
```

## Matrix normalizations

| Function   | Description |
|:---------- |:----------- |
| [`det1`](@ref) | Normalize the determinant|
| [`tr1`](@ref) | Normalize the trace|
| [`normalizeCol!`](@ref) | Normalize one or more columns|

⋅

```@docs
det1
tr1
normalizeCol!
```

## Boolean functions of matrices

| Function   | Description |
|:---------- |:----------- |
| [`ispos`](@ref) | Check whether a real vector or diagonal matrix are comprised of all positive elements|

```@docs
ispos
```

## Scalar functions of matrices

| Function   | Description |
|:---------- |:----------- |
| [`colProd`](@ref) | Sum of products of the elements in two columns |
| [`sumOfSqr`](@ref), `ss` | Sum of squares of all elements or of specified columns |
| [`sumOfSqrDiag`](@ref), `ssd` | Sum of squares of the diagonal elements |
| [`colNorm`](@ref) | Eucliden norm of a column |
| [`sumOfSqrTril`](@ref), `sst` | Sum of squares of the lower triangle elements up to a given underdiagonal |
| [`tr`](@ref) | Fast trace of the product of two Hermitian matrices |
| [`quadraticForm`](@ref), `qf` | Fast quadratic form |
| [`fidelity`](@ref) | (Quantum) Fidelity of two positive matrices |

⋅

```@docs
colProd
sumOfSqr
sumOfSqrDiag
colNorm
sumOfSqrTril
tr
quadraticForm
fidelity
```

## Diagonal functions of matrices

| Function   | Description |
|:---------- |:----------- |
| [`fDiag`](@ref), `𝑓𝔻` | Elemen-wise functions of matrix diagonals|
| [`DiagOfProd`](@ref), `dop` | Diagonal of the product of two matrices|

⋅

```@docs
fDiag
DiagOfProd
```

## Unitary functions of matrices

| Function   | Description |
|:---------- |:----------- |
| [`mgs`](@ref) | Modified Gram-Schmidt orthogonalization|

⋅

```@docs
mgs
```

## Matrix function of matrices

| Function   | Description |
|:---------- |:----------- |
| [`fVec`](@ref) | General function for multi-threaded computation of means and sums of matrix vectors|
| [`congruence`](@ref), `cong` | Compute congruent transformations |




⋅

```@docs
fVec
congruence
```

## Spectral decompositions of positive matrices

| Function   | Description |
|:---------- |:----------- |
| [`evd`](@ref) | Eigenvalue-Eigenvector decomposition of a matrix in ``UΛU'=P`` form|
| [`spectralFunctions`](@ref) | Mother function for creating spectral functions of eigenvalues|
| [`pow`](@ref)| Power of a positive matrix for any number of exponents in one pass|
| [`invsqrt`](@ref)| Principal square root inverse (whitening) of a positive matrix|
| [`sqr`](@ref)| Square of a positive matrix|
| [`powerIterations`](@ref), `powIter` | Power method for estimating any number of eigenvectors and associated eigenvalues|

⋅

```@docs
evd
spectralFunctions
pow
invsqrt
sqr
powerIterations
```

## Decompositions involving triangular matrices

| Function   | Description |
|:---------- |:----------- |
| [`choL`](@ref) | Lower triangula factor of Cholesky decomposition|

⋅

```@docs
choL
```
