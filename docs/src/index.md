# PosDefManifold Documentation

## Requirements

Julia version ≥ 1.0.3

## Installation

Execute the following command in Julia's REPL:

    ]add PosDefManifold

To obtain the latest development version execute instead

    ]add PosDefManifold#master


## Overview

![Figure 1](assets/Fig1.jpg)

**Riemannian geometry** studies smooth manifolds, multi-dimensional curved spaces with peculiar geometries endowed with non-Euclidean metrics. In these spaces Riemannian geometry allows the definition of **angles**, **geodesics** (shortest path between two points), **distances** between points, **centers of mass** of several points, etc.

In this package we are concerned with the manifold **P** of **positive definite matrices**, either *symmetric positive definite* or *Hermitian positive definite*.

In several fields of research such as *computer vision* and *brain-computer interface*, treating data in the **P** manifold has allowed the introduction of machine learning approaches with remarkable characteristics, such as simplicity of use, excellent classification accuracy, as demonstrated by the [winning score](http://alexandre.barachant.org/challenges/) obtained in six international data classification competitions, and the ability to operate transfer learning (Congedo et *al.*, 2017)[🎓](@ref)).

For a formal introduction to the **P** manifold the reader is referred to the monography written by Bhatia (2007)[🎓](@ref).

For an introduction to Riemannian geometry and an overview of mathematical tools implemented in this package, see [Intro to Riemannian Geometry](@ref) in this documentation.

For starting using this package, browse the code units listed here below and execute the many **code examples** you will find therein. The core functions are contained in unit *riemannianGeometry.jl*.

## Code units

**PosDefManifold** includes six code units (.jl files):

| Unit   | Description |
|:----------:| ----------- |
| [MainModule (PosDefManifold.jl)](@ref) | Main module, constants, types, aliases, tips & tricks |
| [riemannianGeometry.jl](@ref) | The fundamental unit collecting all functions acting on the **P** manifold |
| [linearAlgebra.jl](@ref) | Collection of linear algebra routines |
| [signalProcessing.jl](@ref) | Collection of signal processing routines |
| [classification.jl](@ref) | Collection of classification routines |
| [test.jl](@ref) | Unit performing all tests |

## Contents

```@contents
Pages = [       "index.md",
                "introToRiemannianGeometry.md",
                "MainModule.md",
                "riemannianGeometry.md",
                "linearAlgebra.md",
                "signalProcessing.md",
                "classification.md",
                "test.md"]
Depth = 1
```

## Index

```@index
```
