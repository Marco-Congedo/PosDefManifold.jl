# PosDefManifold.jl

| **Documentation**  | 
|:---------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://marco-congedo.github.io/PosDefManifold.jl/dev/) |

**PosDefManifold** is a [**Julia**](https://julialang.org/) package for manipulating data in the manifold **P** of real or complex [**positive definite matrices**](https://en.wikipedia.org/wiki/Definiteness_of_a_matrix). The package supports **10 metrics** acting on **P**, two of which define [**Riemannian manifolds**](https://en.wikipedia.org/wiki/Riemannian_manifold).

![](/docs/src/assets/Fig1.jpg)

This package computes **distances**, **geodesics**, **weighted Fréchet means** and **inductive Fréchet means** for the supported metrics, **graph Laplacians**, **Lapacian eigen maps**, **tangent space projections**, **parallel transport**, the **geometric midrange**, matrix **approximations** and a special solution of the [**Procrustes problem**](https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem) in **P**. Key functions are *multi-threaded*.

The [documentation](https://marco-congedo.github.io/PosDefManifold.jl/dev/) is rich and gives all mathematical details pertaining to the implemented functions.

For **machine learning** in **P** see [PosDefManifoldML.jl](https://github.com/Marco-Congedo/PosDefManifoldML.jl).

For **optimization** in **P** see [Manopt.jl](http://www.manoptjl.org/stable/).

For similar code resources in other programming languages see [here](https://sites.google.com/site/marcocongedo/science/code-resources).

## About the author

[Marco Congedo](https://sites.google.com/site/marcocongedo) is
a research scientist of [CNRS](http://www.cnrs.fr/en) (Centre National de la Recherche Scientifique), working at [UGA](https://www.univ-grenoble-alpes.fr/english/) (University of Grenoble Alpes), in Grenoble (France).

## Contact

marco *dot* congedo *at* gmail *dot* com

## Examples:

```
using PosDefManifold
 
P, Q = randP(20), randP(20)       # random Positive Definite Matrices (PDMs) of size 20x20
d = distance(Fisher, P, Q)        # distance between P and Q 
R = geodesic(Fisher, P, Q, 0.1)   # move on the geodesic relying P and Q 
G = mean(Fisher, P, Q)            # mean of P and Q (geodesic mid-point) 
U = procrustes(P, Q)              # solution to a special Procrustes problem
S = logMap(Fisher, P, G)          # tangent space mapping  
 
𝐏 = randP(20, 100)                 # random set 𝐏 of 100 PDMs of size 20x20
G = mean(Fisher, 𝐏)                # mean of all matrices in set 𝐏 
λ, U, i, c = spEmb(Fisher, 𝐏, 3)   # spectral embedding in 3D
```

