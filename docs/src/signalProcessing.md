# signalProcessing.jl

This unit contains miscellaneous signal processing functions
useful in relation to the Riemannian geometry of the manifold
of *Symmetric Positive Definite (SPD)* or *Hermitian Positive Definite (HPD)*
matrices. In Julia those are `Hermitian` matrices, see [typecasting matrices](@ref).


| Function   | Description |
|:----------:| ----------- |
| [`randChi²`](@ref), randχ²| Generate a random variable distributed as a chi-squared |
| [`randEigvals`](@ref), `randλ` | Generate a random vectors of real positive eigenvalues |
| [`randEigvalsMat`](@ref), `randΛ` | Generate a random diagonal matrix of real positive eigenvalues |
| [`randUnitaryMat`](@ref), `randU` | Generate a random orthogonal or unitary matrix |
| [`randPosDefMat`](@ref), `randP` | Generate one or an array of random positive definite matrices |
| [`regularize!`](@ref) | Regularize an array of positive definite matrices |
| [`gram`](@ref) | Gram matrix of a matrix|
| [`trade`](@ref) | trace and determinant of a matrix as a 2-tuple|
⋅

```@docs
randChi²
randEigvals
randEigvalsMat
randUnitaryMat
randPosDefMat
regularize!
gram
trade
```
