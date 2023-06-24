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
|`sqrt2inv`|1/√2 | 0.7071067811865475 |
|`golden`| (√5+1)/2 | 1.618033988749... |
|`goldeninv`| (√5-1)/2 | 0.618033988749... |
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


All packages above are built-in julia packages.

### types

#### Metric::Enumerated type
```julia
@enum Metric begin
  Euclidean    =1
  invEuclidean =2
  ChoEuclidean =3
  logEuclidean =4
<<<<<<< HEAD
  LogCholesky  =5
=======
  logCholesky  =5
>>>>>>> dev
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
 the metric, for example:

```julia
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

```julia
d=distance(metric, P).
```

To know what is the current metric, you can get it as a string using:

```julia
s=string(metric)
```

To see the list of metrics in type `Metric` use:

```julia
instances(Metric)
```

#### Array of Matrices types

#### 𝕄Vector type
   `𝕄Vector=Vector{𝕄}`

 This is a vector of general `Matrix` matrices, alias of `MatrixVector`.
 Julia sees it as: `Array{Array{T,2} where T,1}`. See [aliases](@ref) for the 𝕄 symbol and [typecasting matrices](@ref) for the use of matrices in **PosDefManifold**.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

 **𝕄Vector₂ type**

   `𝕄Vector₂=Vector{𝕄Vector}`

 This is a vector of [𝕄Vector type](@ref) objects,
 i.e., a vector of vectors of matrices.
 It is the alias of `MatrixVector₂`.
 Julia sees it as: `Array{Array{Array{T,2} where T,1},1}`.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.
    However the several `𝕄Vector` objects it holds do not need to have the same length.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

#### 𝔻Vector type
   `𝔻Vector=Vector{𝔻}`

 This is a vector of `Diagonal` matrices, alias of `DiagonalVector`.
 Julia sees it as: `Array{Diagonal,1}`. See [aliases](@ref) for the 𝔻 symbol and [typecasting matrices](@ref) for the use of Diagonal matrices in **PosDefManifold**.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds  are of the same dimension.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

 **𝔻Vector₂ type**

   `𝔻Vector₂=Vector{𝔻Vector}`

 This is a vector of [𝔻Vector type](@ref) objects,
 i.e., a vector of vectors of `Diagonal` matrices.
 It is the alias of `DiagonalVector₂`.
 Julia sees it as: `Array{Array{Diagonal,1},1}`.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.
    However the several `𝔻Vector` objects it holds do not need to have the same length.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

#### 𝕃Vector type
   `𝕃Vector=Vector{𝕃}`

 This is a vector of `LowerTriangular` matrices, alias of `LowerTriangularVector`.
 Julia sees it as: `Array{LowerTriangular,1}`. See [aliases](@ref) for the 𝕃 symbol and [typecasting matrices](@ref) for the use of LowerTriangular matrices in **PosDefManifold**.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

 **𝕃Vector₂ type**

   `𝕃Vector₂=Vector{𝕃Vector}`

 This is a vector of [𝕃Vector type](@ref) objects, i.e.,
 a vector of vectors of `LowerTriangular` matrices.
 It is the alias of `LowerTriangularVector₂`.
 Julia sees it as: `Array{Array{LowerTriangular,1},1}`.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.
    However the several `𝕃Vector` objects it holds do not need to have the same length.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

#### ℍVector type
   `ℍVector=Vector{ℍ}`

 This is a vector of `Hermitian` matrices, alias of `HermitianVector`.
 Julia sees is at: `Array{Hermitian,1}`.See [aliases](@ref) for the ℍ symbol and [typecasting matrices](@ref) for the use of Hermitian matrices in **PosDefManifold**.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

 **ℍVector₂ type**

    `ℍVector₂=Vector{ℍVector}`

 This is a vector of [ℍVector type](@ref)
 objects, i.e., a vector of vectors of `Hermitian` matrices.
 It is the alias of `HermitianVector₂`.
 Julia sees it as: `Array{Array{Hermitian,1},1}`.

!!! warning "Nota bene"
    This object is meant to hold matrices living in the same manifold,
    therefore it is assumed by all methods that all matrices it holds are of the same dimension.
    However the several `ℍVector` objects it holds do not need to have the same length.

 **See** [`dim`](@ref), [`typeofMatrix`](@ref)

#### RealOrComplex type
   `RealOrComplex=Union{Real, Complex}`

 This is the [Union](https://docs.julialang.org/en/v1/base/base/#Core.Union) of `Real` and `Complex` types.

#### AnyMatrix type
   `AnyMatrix=Union{𝔻{T}, 𝕃{T}, ℍ{T}, 𝕄{T}} where T<:RealOrComplex`

 This is the [Union](https://docs.julialang.org/en/v1/base/base/#Core.Union)
 of real or complex `Diagonal`, `LowerTriangular`, `Hermitian` and `Matrix` types. It is often used in the definition of functions.

 **See** [aliases](@ref)

#### AnyMatrixVector type
   `AnyMatrixVector=Union{𝕄Vector, 𝔻Vector, 𝕃Vector, ℍVector}`

 This is the [Union](https://docs.julialang.org/en/v1/base/base/#Core.Union) of 𝕄Vector, 𝔻Vector, 𝕃Vector and ℍVector. It is often used in the definition of functions.
 See [Array of Matrices types](@ref).

 **AnyMatrixVector₂ type**

   `AnyMatrixVector₂=Union{𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂, ℍVector₂}`

 This is the [Union](https://docs.julialang.org/en/v1/base/base/#Core.Union) of 𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂, ℍVector₂. It is often used in the definition of functions. See [Array of Matrices types](@ref).

### tips & tricks

#### typecasting matrices
 Several functions in **PosDefManifold** implement multiple dispatch and can handle several kinds of matrices as input, however the core functions for manipulating objects on the Riemannian manifold of positive definite matrices act by definition on positive definite matrices only.
 Those matrices must therefore be either
 *symmetric positive definite (SPD, real)* or *Hermitian positive definite (HPD, complex)*.
 Such matrices are uniformly identified in **PosDefManifold** as being of the `Hermitian` type, using the standard [LinearAlgebra](https://bit.ly/2JOiROX) package.
 The alias `ℍ` is used consistently in the code (see [aliases](@ref)).
 If the input is not flagged as `Hermitian`, the functions restricting the input to *positive definite matrices* will not be accessible.

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

 Be careful: `Hermitian(P)` will construct and Hermitian matrix from the argument.
 If the matrix argument is not symmetric (if real) or Hermitian (if complex)
 it will be made so by copying the transpose (if real) or complex conjugate
 and transpose (if complex) of a triangular part into the other.
 See [Hermitian](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.Hermitian).

 If you want to construct an [ℍVector type](@ref) from, say, two Hermitian matrices `P` and `Q`, don't write `A=[P, Q]`, but rather `A=ℍVector([P, Q])`. In fact,
 the first is seen by Julia as

    2-element Array{Hermitian{Float64,Array{Float64,2}},1},

 while the latter as

    2-element Array{Hermitian,1},

 which is the type expected in all functions taking an `ℍVector type` as argument.

 Other functions act on generic matrices (of type [Matrix](https://docs.julialang.org/en/v1/base/arrays/#Base.Matrix)). This is seen by Julia as `Array{T,2} where T`.
 Keep in mind that the functions writing on the argument matrix such as
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

```julia
Matrix(X')
```

 - here is how to get an `Hermitian` matrix out of the diagonal part of an `Hermitian` matrix H:

```julia
Hermitian(Matrix(Diagonal(H)))
```

 - here is how to get a `LowerTriangular` matrix out of an `Hermitian` matrix H:

```julia
LowerTriangular(Matrix(H))
```

 For example, you can use this to pass a full inter-distance matrix to the [`laplacian`](@ref) function to obtain the Laplacian matrix.

 A useful function is [`typeofMatrix`](@ref). For example, the following line
 typecasts matrix `M` to the type of matrix `P` and put the result in `A`:

```julia
A=typeofMatrix(P)(M)
```

#### Threads
Some functions in **PosDefManifold** explicitly call BLAS routines
for optimal performnce. This is reported in the help section of the
concerned functions. Most functions calls BLAS routine implicitly
via Julia. You can set the number of threads
the BLAS library should use by:

```julia
using LinearAlgebra
BLAS.set_num_threads(n)
```

where `n` is the number of threads.
By default, **PosDefManifold** reserves to BLAS
all CPU threads available on your computer (given by the output of `Sys.CPU_THREADS`).
The number of threads used by Julia
for multi-threaded computations is given by the output of function `Threads.nthreads()`.
In Windows this latter number of threads is set to half the available threads.
In Linux and OSX defaults to one and is controlled by an environment variable, i.e.,

   `export JULIA_NUM_THREADS=4`.

In Linux, working with the Atom IDE, you also have to
set to `global` the field found in Atom under
`Settings(or Preferences)/julia-client/Settings/Julia Options/Number of Threads`.

In Windows, set the desired number of threads in the settings
of the julia-client Juno package.

See for example this [post](https://discourse.julialang.org/t/issue-number-of-threads/14593), this [post](https://discourse.julialang.org/t/customize-number-of-threads-interactively/11574/2) and julia
[doc on threads](https://docs.julialang.org/en/v1/manual/parallel-computing/#Multi-Threading-(Experimental)-1).

Notice that **PosDefManifold** features many multi-threaded functions and these
may allow a gain in computation time only if Julia is instructed to use
at least two threads.
