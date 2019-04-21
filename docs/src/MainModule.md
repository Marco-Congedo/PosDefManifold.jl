# MainModule (PosDefManifold.jl)

This is the main unit containing the **PosDefManifold** *module*.

It uses the following standard Julia packages:

| using  |
|:----------:|
| [Linear Algebra](https://bit.ly/2W5Wq8W) |
| [Statistics](https://bit.ly/2Oem3li) |

Examples in some units of **PosDefManifold** also uses the `Plots` package.

The main module does not contains functions, but it declares all **constant**,
**types** and **aliases** of Julia functions and types used in all units.

| Contents  |
|:----------:|
| [constants](@ref) |
| [aliases](@ref) |
| [types](@ref) |
| [tips & tricks](@ref) |

## constants
| constant   | value | numeric value
|:----------:| ----------- | ----------- |
|`sqrt2` |√2 | 1.4142135623730951 |
|`invsqrt2`|1/√2 | 0.7071067811865475 |
|`minpos`| 1e-15 | 0.000000000000001|
|`maxpos`| 1e15 | 100000000000000|

## aliases

| alias   | Julia function | in Package | tab-completition | REPL support |
|:----------:| ----------- | ----------- | ----------- | ----------- |
|`𝚺` |[`sum`](https://bit.ly/2FcsAJg)|Base| \bfSigma | ⛔ |
|`𝛍`|[`mean`](https://bit.ly/2TOakA0)|Statistics| \bfmu | ⛔ |
|`⋱`|[`Diagonal`](https://bit.ly/2Jovxf8)|LinearAlgebra| \ddots | ✓ |
|`ℍ`|[`Hermitian`](https://bit.ly/2JOiROX)|LinearAlgebra| \bbH | ✓ |

All packages above are built-in julia packages.

## types

### Metric::Enumerated type
```
@enum Metric begin
  Euclidean    =1
  invEuclidean =2
  ChoEuclidean =3
  logEuclidean =4
  LogCholesky  =5
  Fisher       =6 # default metric
  logdet0      =7
  Jeffrey      =8
  VonNeumann   =9
  Wasserstein  =10
end
```

 Riemannian manipulations are defined for a given *metric* (see [metrics](@ref)).
 An instance for this type is requested as an argument in many functions
 contained in the [riemannianGeometry.jl](@ref) unit in order to specify
 the metric, unless the default metric (Fisher) is sought.

```
 ## Example
 # generate a 15x15 symmetric positive definite matrix
 P=randP(15)              
 # distance from P to the identity matrix according to the logdet0 metric
 d=distance(P, logdet0)  
```

 If you want to work consistently with a specific metric,
 you may want to declare in your script a global variable such as

    global metric=logdet0  or  global metric=Metric(Int(logdet0)),

 and then pass `metric` as argument in all your computations,
 *e.g.*, referring to the above example,

    d=distance(P, metric).

To know what is the current metric, get it as a string as:

    s=string(metric)

### RealOrComplex type
 `RealOrComplex=Union{Real, Complex}` is the Union of Real and Complex Types.

### ℍVector type
 `ℍVector=Vector{ℍ}` is a vector of Hermitian matrices.
 Julia sees is at: `Array{Hermitian,1}`.

**ℍVector₂ type**
  `ℍVector₂=Vector{ℍVector}` is a vector of [ℍVector type](@ref) objects, i.e., a vector of vectors of Hermitian matrices. Julia sees it as: `Array{Array{Hermitian,1},1}`. Note that `ℍVector₂` is not a matrix
  of Hermitian matrices since the several `ℍVector` it holds do not need
  to have the same length. See [aliases](@ref) for the ℍ symbol and
  [typecasting matrices](@ref) for the use of Hermitian matrices  in **PosDefManifold**.

## tips & tricks

### typecasting matrices
 Several functions in **PosDefManifold** implement multiple dispatch and can handle  
 several kinds of matrices as input, however the core functions for manipulating  
 objects on the Riemannian manifold of positive definite matrices act by definition
 on positive definite matrices only.
 Those matrices must therefore be either
 *symmetric positive definite (real)* or *Hermitian (complex)*.
 Such matrices are uniformly identified in **PosDefManifold** as being of the `Hermitian` type, using the standard [LinearAlgebra](https://bit.ly/2JOiROX) package.
 The alias `ℍ` is used consistently in the code (see [aliases](@ref)).
 If the input is not flagged, the functions restricting the input to *positive definite matrices* will give an error.

 **Example**

    julia> using LinearAlgebra

    julia> f(S::Hermitian)=S*S'
    f (generic function with 1 method)

    julia> A=randn(3, 3)
    3×3 Array{Float64,2}:
     -0.67407  -0.344258    0.203714
     -1.06551  -0.0233796   0.975465
     -1.04727  -1.19807    -0.0219121

    julia> H=A*A' # although SPD, H is not automatically flagged as Hermitian
    3×3 Array{Float64,2}:
     0.614384  0.924991  1.11391
     0.924991  2.08738   1.12251
     1.11391   1.12251   2.53263

    julia> f(H)
    ERROR: MethodError: no method matching f(::Array{Float64,2})
    Closest candidates are:
      f(::Hermitian) at none:1

 If you construct a positive definite matrix and it is not flagged,
 you can do so simply by **typecasting** it, that is, passing as argument to the
 functions `Hermitian(P)` instead of just `P`. The `ℍ` alias can be
 used for short, *i.e.*, `ℍ(P)`. Continuing the example above:

    julia> f(ℍ(H))  # this way it works, equivalent to f(Hermitian(H))
    3×3 Array{Float64,2}:
     2.47388  3.74948  4.54381
     3.74948  6.4728   6.21635
     4.54381  6.21635  8.91504

 Similarly, if you want to construct an [ℍVector type] from, say, two Hermitian
 matrices `P` and `Q`, don't write `A=[P, Q]`, but rather `A=ℍVector([P, Q])`.

 Other functions act on generic matrices (of type [Matrix](https://docs.julialang.org/en/v1/base/arrays/#Base.Matrix)).
 To those functions you can pass any matrix.
 However, keep in mind that the functions writing on the argument matrix such as
 [`normalizeCol!`](@ref) will give an error if you pass an `Hermitian` matrix,
 since Julia does not allow writing on non-diagonal elements of those matrices.
 In this case typecast it in another object using the `Matrix` type;
 suppose `H` is `Hermitian`, you would use for example:

    julia> X=Matrix(H)
    julia> normalizeCol!(X, 1)
    julia> norm(X[:, 1])
    1.0

 Another example when typecasting is useful: the [`gram`](@ref) function
 takes a `Matrix` type as argument (since `X` is expected to be a data matrix),
 like in

    H=gram(X)

 The following will not work though:
```
 H=gram(X')
```

 since `X'` is an `Adjoint` type. The problem is fixed by typecasting the
 adjoint matrix, such as

  H=gram(Matrix(X'))

 Another example: here is how to get an Hermitian matrix out of the
 diagonal part of an Hermitian matrix H:

    Hermitian(Matrix(Diagonal(H))).
