# MainModule (PosDefManifold.jl)

This is the main unit containing the **PosDefManifold** *module*.

It uses the following standard Julia packages:

| using  |
|:----------:|
| [LinearAlgebra](https://bit.ly/2W5Wq8W) |
| [Statistics](https://bit.ly/2Oem3li) |

Examples in some units of **PosDefManifold** also uses the `Plots` package. Take a look at this [tutorial](https://cleytonfar.github.io/posts/using-julia-for-data-science-part-03/)
for an introduction to data plotting with Julia.

The main module does not contains functions, but it declares all **constant**,
**types** and **aliases** of Julia functions used in all units.

| Contents  |
|:----------:|
| [constants](@ref) |
| [aliases](@ref) |
| [types](@ref) |
| [tips & tricks](@ref) |

### constants
| constant   | value | numeric value
|:----------:| ----------- | ----------- |
|`sqrt2` |√2 | 1.4142135623730951 |
|`invsqrt2`|1/√2 | 0.7071067811865475 |
|`maxpos`| 1e15 | 100000000000000|


### aliases

| alias   | Julia function | in Package | tab-completition | REPL support |
|:----------:| ----------- | ----------- | ----------- | ----------- |
|`𝚺` |[`sum`](https://bit.ly/2FcsAJg)|Base| \bfSigma | ⛔ |
|`𝛍`|[`mean`](https://bit.ly/2TOakA0)|Statistics| \bfmu | ⛔ |
|`𝕄`|[`Matrix`](https://docs.julialang.org/en/v1/base/arrays/#Base.Matrix)|Base| \bbM | ⛔ |
|`𝔻`|[`Diagonal`](https://bit.ly/2Jovxf8)|LinearAlgebra| \bbD | ⛔ |
|`ℍ`|[`Hermitian`](https://bit.ly/2JOiROX)|LinearAlgebra| \bbH | ✓ |
|`𝕃`|[`LowerTriangular`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.LowerTriangular)|LinearAlgebra| \bbH | ⛔ |
|`❗`|[`Threads.@threads`](https://docs.julialang.org/en/v1/base/multi-threading/#Base.Threads.@threads)|Base.Threads| | ✓ |



All packages above are built-in julia packages.

### types

#### Metric::Enumerated type
```
@enum Metric begin
  Euclidean    =1
  invEuclidean =2
  ChoEuclidean =3
  logEuclidean =4
  LogCholesky  =5
  Fisher       =6
  logdet0      =7
  Jeffrey      =8
  VonNeumann   =9
  Wasserstein  =10
end
```

 Riemannian manipulations are defined for a given *metric* (see [metrics](@ref)).
 An instance for this type is requested as an argument in many functions
 contained in the [riemannianGeometry.jl](@ref) unit in order to specify
 the metric.

```
 ## Example
 # generate a 15x15 symmetric positive definite matrix
 P=randP(15)              
 # distance from P to the identity matrix according to the logdet0 metric
 d=distance(logdet0, P)  
```

 If you want to work consistently with a specific metric,
 you may want to declare in your script a global variable such as

    global metric=logdet0  or  global metric=Metric(Int(logdet0)),

 and then pass `metric` as argument in all your computations,
 *e.g.*, referring to the above example,

    d=distance(metric, P).

To know what is the current metric, you can get it as a string using:

    s=string(metric)

To see the list of metrics in type Metric use:

    instances(Metric)

#### RealOrComplex type
 `RealOrComplex=Union{Real, Complex}`

 This is the Union of Real and Complex Types.

#### 𝔻Vector type
  `𝔻Vector=Vector{𝔻}`

  This is a vector of `Diagonal` matrices, alias of `DiagonalVector`.
  Julia sees is at: `Array{Diagonal,1}`.See [aliases](@ref) for the 𝔻 symbol and [typecasting matrices](@ref) for the use of Diagonal matrices in **PosDefManifold**.

 **𝔻Vector₂ type**

 `𝔻Vector₂=Vector{𝔻Vector}` is a vector of [𝔻Vector type](@ref) objects, i.e., a vector of vectors of `Diagonal` matrices.
 It is the alias of `DiagonalVector₂`.
 Julia sees it as: `Array{Array{Diagonal,1},1}`. Note that `𝔻Vector₂` is not a matrix of Hermitian matrices since the several `ℍVector` objects it holds do not need to have the same length.

#### ℍVector type
 `ℍVector=Vector{ℍ}`

 This is a vector of `Hermitian` matrices, alias of `HermitianVector`.
 Julia sees is at: `Array{Hermitian,1}`.See [aliases](@ref) for the ℍ symbol and [typecasting matrices](@ref) for the use of Hermitian matrices in **PosDefManifold**.

**ℍVector₂ type**

 `ℍVector₂=Vector{ℍVector}` is a vector of [ℍVector type](@ref) objects, i.e., a vector of vectors of `Hermitian` matrices.
 It is the alias of `HermitianVector₂`.
 Julia sees it as: `Array{Array{Hermitian,1},1}`. Note that `ℍVector₂` is not a matrix of Hermitian matrices since the several `ℍVector` objects it holds do not need to have the same length.

### tips & tricks

#### typecasting matrices
 Several functions in **PosDefManifold** implement multiple dispatch and can handle several kinds of matrices as input, however the core functions for manipulating objects on the Riemannian manifold of positive definite matrices act by definition on positive definite matrices only.
 Those matrices must therefore be either
 *symmetric positive definite (SPD, real)* or *Hermitian positive definite (HPD, complex)*.
 Such matrices are uniformly identified in **PosDefManifold** as being of the `Hermitian` type, using the standard [LinearAlgebra](https://bit.ly/2JOiROX) package.
 The alias `ℍ` is used consistently in the code (see [aliases](@ref)).
 If the input is not flagged as `Hermitian`, the functions restricting the input to *positive definite matrices* will give an error.

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

 Similarly, if you want to construct an [ℍVector type](@ref) from, say, two Hermitian matrices `P` and `Q`, don't write `A=[P, Q]`, but rather `A=ℍVector([P, Q])`. In fact,
 the first is seen by Julia as

    2-element Array{Hermitian{Float64,Array{Float64,2}},1},

 while the latter as

    2-element Array{Hermitian,1},

 which is the type expected in all functions taking an `ℍVector type` as argument.

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

  Some more examples:


 - Typecasting `Adjoint` matrices:

    Matrix(X')

 - here is how to get an `Hermitian` matrix out of the
 diagonal part of an `Hermitian` matrix H:

    Hermitian(Matrix(Diagonal(H))).

 - here is how to get a `LowerTriangular` matrix out of an
 `Hermitian` matrix H:

    LowerTriangular(Matrix(H)).

 For example, you can use this to pass a full inter-distance matrix to the [`laplacian`](@ref) function to obtain the Laplacian matrix.


#### BLAS routines
Some functions in **PosDefManifold** call BLAS routines for optimal performnce.
This is reported in the help section of the concerned functions.
In most functions julia calls BLAS routines automatically.
You can set the number of threads
the BLAS library should use by:

    using LinearAlgebra
    BLAS.set_num_threads(n)

where `n` is the number of threads.
By default, **PosDefManifold** reserves to BLAS
all CPU threads available on your computer (i.e., the output of `Sys.CPU_THREADS`) minus Threads.nthreads().
See this [post](https://discourse.julialang.org/t/issue-number-of-threads/14593), this [post](https://discourse.julialang.org/t/customize-number-of-threads-interactively/11574/2) and julia
[doc on threads](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1).
