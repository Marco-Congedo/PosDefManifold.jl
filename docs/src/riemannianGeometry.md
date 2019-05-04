# riemannianGeometry.jl

This is the fundamental unit of **PosDefManifold**. It contains functions
for manipulating points in the Riemannian manifold of
*Symmetric Positive Definite (SPD)* or *Hermitian Positive Definite (HPD)* matrices. In Julia those are `Hermitian` matrices, see [typecasting matrices](@ref).

The functions are divided in six categories:

| Category  | Output |
|:----------:| ----------- |
| 1. [Geodesic equations](@ref) | interpolation, extrapolation,... |
| 2. [Distances](@ref) | length of geodesics |
| 3. [Graphs and Laplacians](@ref) | for spectral embedding, eigenmaps, system dynamics,...|
| 4. [Means](@ref) | mid-points of geodesics, centers of mass of several points |
| 5. [Tangent Space operations](@ref) | maps from the manifold to the tangent space and viceversa |
| 6. [Procrustes problems](@ref) | for data matching, transfer learning,...|

⋅

## Geodesic equations

```@docs
geodesic
```

## Distances

```@docs
distanceSqr
distance
```

## Graphs and Laplacians
```@docs
distanceSqrMat
distanceMat
laplacian
laplacianEigenMaps
spectralEmbedding
```

## Means
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
```@docs
logMap
expMap
vecP
matP
```

## Procrustes problems
```@docs
procrustes
```
