# riemannianGeometry.jl

This is the fundamental unit of **PosDefManifold**. It contains functions
for manipulating points in the Riemannian manifold of
*Symmetric Positive Definite (SPD)* or *Hermitian Positive Definite (HPD)* matrices. In Julia those are `Hermitian` matrices, see [typecasting matrices](@ref).

The functions are divided in six categories:

| Category  | Output |
|:---------- |:----------- |
| 1. [Geodesic equations](@ref) | interpolation, extrapolation,... |
| 2. [Distances](@ref) | length of geodesics |
| 3. [Graphs and Laplacians](@ref) | for spectral embedding, eigenmaps, system dynamics,...|
| 4. [Means](@ref) | mid-points of geodesics, centers of mass of several points |
| 5. [Tangent Space operations](@ref) | maps from the manifold to the tangent space and viceversa |
| 6. [Procrustes problems](@ref) | for data matching, transfer learning,...|

⋅

## Geodesic equations

| Function   | Description |
|:---------- |:----------- |
| [`geodesic`](@ref) | Geodesic equations (weighted mean of two positive definite matrices) for any metric |

⋅

```@docs
geodesic
```

## Distances

| Function   | Description |
|:---------- |:----------- |
| [`distanceSqr`](@ref), `distance²` | Squared distance between positive definite matrices|
| [`distance`](@ref) | Distance between positive definite matrices|

⋅

```@docs
distanceSqr
distance
```

## Graphs and Laplacians

| Function   | Description |
|:---------- |:----------- |
| [`distanceSqrMat`](@ref), `distance²Mat` | Lower triangular matrix of all squared inter-distances|
| [`distanceMat`](@ref) | Lower triangular matrix of all inter-distances|
| [`laplacian`](@ref) | Laplacian of a squared inter-distances matrix|
| [`laplacianEigenMaps`](@ref), `laplacianEM` | Eigen maps (eigenvectors) of a Laplacian|
| [`spectralEmbedding`](@ref) | Spectral Embedding (all the above function run in series)|

⋅

```@docs
distanceSqrMat
distanceMat
laplacian
laplacianEigenMaps
spectralEmbedding
```

## Means

| Function   | Description |
|:---------- |:----------- |
| [`mean`](@ref) | Weighted Fréchet mean (wFm) of a matrix set using any metric |
| [`means`](@ref) | As above for several sets at once |
| [`generalizedMean`](@ref) | Generalized wFm of a matrix set |
| [`geometricMean`](@ref), `gMean` | wFm of a matrix set according to the Fisher metric (iterative)|
| [`logdet0Mean`](@ref), `ld0Mean` | wFm of a matrix set according to the logdet0 metric (iterative)|
| [`wasMean`](@ref) | wFm of a matrix set according to the Wasserstein metric (iterative)|
| [`powerMean`](@ref) | Power wFm of a matrix set (iterative)|

⋅

```@docs
mean
means
generalizedMean
geometricMean
logdet0Mean
wasMean
powerMean
```

## Tangent Space operations

| Function   | Description |
|:---------- |:----------- |
| [`logMap`](@ref) | Logarithmic map (from manifold to tangent space) |
| [`expMap`](@ref) | Exponential map (from tangent space to manifold) |
| [`vecP`](@ref) | vectorization of matrices in the tangent space |
| [`matP`](@ref) | matrization of matrices in the tangent space (inverse of ``vecp`)|

⋅

```@docs
logMap
expMap
vecP
matP
```

## Procrustes problems

| Function   | Description |
|:---------- |:----------- |
| [`procrustes`](@ref) | Solution to the Procrustes problem in the manifold of positive definite matrices |

⋅

```@docs
procrustes
```
