# PosDefManifold Documentation

## Requirement

Julia v1.0.3 or higher

## Installation

The package is not registered yet. Use the [github
url](https://github.com/Marco-Congedo/PosDefManifold.jl/tree/master/src) for
adding the **PosDefManifold** package.

## Overview

![Figure 1](assets/Fig1.jpg)

**Riemannian geometry** studies smooth manifolds, multi-dimensional curved spaces with peculiar geometries endowed with non-Euclidean metrics. In these spaces Riemannian geometry allows the definition of **angles**, **geodesics** (shortest path between two points), **distances** between points, **centers of mass** of several points, etc.

In this package we are concerned with the manifold **P** of **positive definite matrices**, either *symmetric positive definite* or *Hermitian positive definite*.

In several fields of research such as *computer vision* and *brain-computer interface*, treating data in the **P** manifold has allowed the introduction of machine learning approaches with remarkable characteristics, such us simplicity of use, excellent classification accuracy, as demonstrated by the [winning score](http://alexandre.barachant.org/challenges/) obtained in six international data classification competitions, and the ability to operate transfer learning (Congedo et *al.*, 2017)[🎓](@ref)).

For a formal introduction to the **P** manifold the reader is referred to Bhatia (2007)[🎓](@ref).

For an introduction to Riemannian geometry and an overview of mathematical tools implemented in this package, see [Intro To Riemannian Geometry](introToRiemannianGeometry/index.html) in this documentation.

For starting using this package, browse the code units listed here below and execute the many **example code** you will find therein.

## Code units

**PosDefManifold** includes five code units (.jl files):

| Unit   | Description |
|:----------:| ----------- |
| [PosDefManifold.jl](MainModule/index.html) | Main module |
| [riemannianGeometry.jl](riemannianGeometry/index.html) | The fundamental unit collecting all functions acting on the **P** manifold |
| [linearAlgebra.jl](linearAlgebra/index.html) | Collection of linear algebra routines |
| [signalProcessing.jl](signalProcessing/index.html) | Collection of signal processing routines |
| [test.jl](test/index.html) | Unit performing all tests |

## Contents

```@contents
Pages = [ "index.md", "introToRiemannianGeometry.md", "MainModule.md", "riemannianGeometry.md", "linearAlgebra.md", "signalProcessing.md", "test.md"]
Depth = 1
```

## Index

```@index
```
