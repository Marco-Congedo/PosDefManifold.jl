#    Unit signalProcessing.jl, part of PosDefManifold Package for julia language
#
#    MIT License
#    Copyright (c) 2019-22, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home
#
#    DESCRIPTION
#    This Unit implements Signal Processing methods and prcedures
#    useful in relation to the Riemannian Geometry of the manifold
#    of Symmetric Positive Definite (SPD) or Hermitian matrices


"""
    randChi²(df::Int)

 **alias**: `randχ²`

 Generate a random variable distributed as a *chi-squared* with `df`
 degrees of freedom.

 It uses the *Wilson–Hilferty transformation* for `df`>=20 -
 see [chi-squared distribution](https://en.wikipedia.org/wiki/Chi-squared_distribution).

 **Examples**
```julia
using Plots, PosDefManifold
chi=[randχ²(2) for i=1:10000]
histogram(chi) # needs Plots package. Check your plots back-end.
```

"""
randChi²(df::Int) =
    df<20 ? sum(randn()^2 for i=1:df) : df*(1.0-2.0/(9.0*df)+randn()*sqrt2/sqrt(9.0*df))^3

randχ²=randChi²

"""
```
    randEigvals(n::Int;
    <
    df::Int=2,
    eigvalsSNR::Real=10e3 >)
```

 **alias**: `randλ`

 Generate an ``n``-vector of random real positive eigenvalues.
 The eigenvalues are generated as in function `randΛ`([`randEigvalsMat`](@ref)),
 the syntax of which is used.

 **See also**: `randU` ([`randUnitaryMat`](@ref)), `randP` ([`randPosDefMat`](@ref)).

 **Examples**
```julia
using Plots, PosDefManifold
λ=sort(randλ(10), rev=true)
σ=sort(randλ(10, eigvalsSNR=10), rev=true)
plot(λ) # needs Plots package. Check your plots back-end.
plot!(σ) # needs Plots package. Check your plots back-end.
```

"""
randEigvals(n::Int;
            df::Int=2,
            eigvalsSNR::Real=10e3) =
    eigvalsSNR==Inf ? [randχ²(df) for i in 1:n] : [randχ²(df)+(df/eigvalsSNR) for i in 1:n]

randλ=randEigvals

"""
```
    (1) randEigvalsMat(n::Int;
    <
    df::Int=2,
    eigvalsSNR::Real=10e3 >)

    (2) randEigvalsMat(n::Int, k::Int;
    < same keyword arguments as in (1) >)
```

 **alias**: `randΛ`

 (1) Generate an ``n⋅n`` diagonal matrix of random real positive eigenvalues.

 (2) An array 1d (of [𝔻Vector type](@ref)) of ``k`` matrices of the kind in (1)

 The eigenvalues are generated according to model

 ``λ_i=χ_{df}^2+η,\\hspace{6pt}\\textrm{for}\\hspace{2pt}i=1:n,``

 where
 - ``χ_{df}^2`` (signal term) is randomly distributed as a [chi-square](https://bit.ly/1IXkulE) with `df` degrees of freedom,
 - ``η`` is a [white noise](https://bit.ly/2TN8472) term, function of *<keyword argument>* `eigvalsSNR`, such that

 ``\\textrm{eigenvalues SNR}=\\mathbb{E}\\big(\\sum_{i=1}^{n}λ_i\\big)\\big/nη.``

 The expected sum ``\\mathbb{E}\\big(\\sum_{i=1}^{n}λ_i\\big)`` here above is the
 expected variance of the signal term, i.e., ``n(df)``, since the expectation
 of a random chi-squared variable is equal to its degrees of freedom.

 If `eigvalsSNR=Inf` is passed as argument, then ``η`` is set to zero, *i.e.*,
 no white noise is added. In any case `eigvalsSNR` must be positive.

 Note that with the default value of *<keyword argument>* `df` (`df=2`)
 the generating model assumes that the eigenvalues
 have exponentially decaying variance, which is often observed on real data.

!!! note "Nota Bene"
    The *<keyword argument>* `eigvalsSNR` expresses the expected eigenvalues SNR
    ([signal-to-noise ratio](https://bit.ly/1VvpvnQ)),
    not the real one, and is not expressed in decibels,
    but as the expected SNR variance ratio.

 This function is used by function `randP` ([`randPosDefMat`](@ref)) to generate
 random positive definite matrices with added white noise in order
 to emulate eigenvalues observed in real data and to
 improve the conditioning of the generated matrices with respect to inversion.



 **See also**: `randλ` ([`randEigvals`](@ref)), `randU` ([`randUnitaryMat`](@ref)),
 `randP` ([`randPosDefMat`](@ref)), `randχ²` ([`randChi²`](@ref)).

 **Examples**
```julia
using PosDefManifold
# (1)
n=3;
U=randU(n);
Λ=randΛ(n, eigvalsSNR=100)
P=U*Λ*U' # generate an SPD matrix
using LinearAlgebra
Q=ℍ(U*Λ*U') # generate an SPD matrix and flag it as 'Hermitian'

# (2) generate an array of 10 matrices of simulated eigenvalues
Dvec=randΛ(n, 10)
```

"""
randEigvalsMat(n::Int;
                df::Int=2,
                eigvalsSNR::Real=10e3) =
        𝔻(randλ(n, df=df, eigvalsSNR=eigvalsSNR))

randEigvalsMat(n::Int, k::Int;
                df::Int=2,
                eigvalsSNR::Real=10e3) =
        𝔻Vector([randΛ(n, df=df, eigvalsSNR=eigvalsSNR) for j=1:k])

randΛ=randEigvalsMat


"""
    (1) randUnitaryMat(n::Int)
    (2) randUnitaryMat(::Type{Complex{T}}, n::Int)

 **aliases**: `randOrthMat`, `randU`

 Generate a random ``n⋅n``
 - (1) [orthogonal](https://bit.ly/2vrr0wU) matrix (real)
 - (2) [unitary](https://bit.ly/2JCHbmC) matrix (complex)

 The matrices are generated running the modified (stabilized)
 [Gram-Schmidt orthogonalization](https://bit.ly/2YE6zvy)
 procedure ([`mgs`](@ref)) on an ``n⋅n`` matrix filled with random Gaussian elements.

 **See also**: `randΛ` ([`randEigvals`](@ref)), `randP` ([`randPosDefMat`](@ref)).

 **Examples**
```julia
using PosDefManifold
n=3;
X=randU(n)*sqrt(randΛ(n))*randU(n)'  # (1) generate a random square real matrix

U=randU(ComplexF64, n);
V=randU(ComplexF64, n);
Y=U*sqrt(randΛ(n))*V' # (2) generate a random square complex matrix
```

"""
randUnitaryMat(n::Int)=mgs(randn(Float64, n, n))

randOrthMat(n::Int)=mgs(randn(Float64, n, n))

randUnitaryMat(::Type{Complex{T}}, n::Int) where {T<:AbstractFloat} = mgs(randn(ComplexF64, n, n))

randU=randUnitaryMat


"""
```
    (1) randPosDefMat(n::Int;
    <
    df::Int=2,
    eigvalsSNR::Real=10e3 >)

    (2) randPosDefMat(::Type{Complex{T}}, n:: Int;
    < same keyword arguments as in (1) >)

    (3) randPosDefMat(n::Int, k::Int;
    <
    df::Int=2,
    eigvalsSNR::Real=10e3,
    SNR::Real=100,
    commuting=false >)

    (4) randPosDefMat(::Type{Complex{T}}, n::Int, k::Int;
    < same keyword arguments as in (3) >)
```

 **alias**: `randP`

 Generate
 - (1) one random `Hermitian` positive definite matrix (real) of size ``n⋅n``
 - (2) one random `Hermitian` positive definite matrix (complex) of size ``n⋅n``
 - (3) an array 1d (of [ℍVector type](@ref)) of ``k`` matrices of the kind in (1)
 - (4) an array 1d (of [ℍVector type](@ref)) of ``k`` matrices of the kind in (2).

 Methods (3) and (4) are **multi-threaded**. See [Threads](@ref).

 For (1) and (2) the matrix is generated according to model

 ``UΛU^H+ηI``,

 where ``U`` is a random orthogonal (1) or unitary (2) matrix generated by
 function `randU`([`randUnitaryMat`](@ref)) and ``Λ``, ``η`` are a positive definite
 diagonal matrix and a non-negative scalar depending on *<optional keywords arguments>*
 `df` and `eigvalsSNR` randomly generated calling function
 `randΛ`([`randEigvalsMat`](@ref)).

 For (3) and (4), if the *<optional keyword argument>* `commuting=true` is passed,
 the ``k`` matrices are generated according to model

  ``UΛ_iU^H+ηI,\\hspace{8pt}``, for ``i``=1:``k``

 otherwise they are generated according to model

 ``(UΛ_iU^H+ηI)+φ(V_iΔ_iV_i^H+ηI),\\hspace{8pt}``, for ``i``=1:``k``  Eq.[1]

 where
 - ``U`` and the ``V_i`` are random (3) orthogonal/(4) unitary matrices,
 - ``Λ_i`` and ``Δ_i`` are positive definite diagonal matrices
 - ``η`` is a non-negative scalar.
  All variables here above are randomly generated as in (1) and (2)
 - ``φ`` is adjusted so as to obtain a desired output `SNR` ([signal-to-noise ratio](https://bit.ly/1VvpvnQ)), which is also an
 *<optional keywords arguments>*, such as

 ``SNR=\\frac{\\displaystyle\\sum_{i=1}^{k}\\textrm{tr}(UΛ_iU^H+ηI)}{\\displaystyle\\sum_{i=1}^{k}\\textrm{tr}φ(V_iΔ_iV_i^H+ηI)}``.

!!! note "Nota Bene"
    The keyword arguments `SNR` is not expressed in decibels,
    but as the expected SNR variance ratio. It must be a positive number.

 A slightly different version of this model for generating positive definite
 matrices has been proposed in (Congedo et *al.*, 2017b)[🎓];
 in the model of Eq. [1]
 - ``UΛ_iU^H`` is the signal term, where the signal is supposed sharing the same coordinates for all matrices,
 - ``φ(V_iΔ_iV_i^H)`` is a structured noise term, which is different for all matrices
 - ``ηI`` is a [white noise](https://bit.ly/2TN8472) term, with same variance for all matrices.

 **See also**: the aforementioned paper and `randΛ` ([`randEigvalsMat`](@ref)).

 **Examples**
```julia
using PosDefManifold
R=randP(10, df=10, eigvalsSNR=1000) # 1 SDP Matrix of size 10x10 #(1)
H=randP(ComplexF64, 5, eigvalsSNR=10) # 1 Hermitian Matrix of size 5x5 # (2)
ℛ=randP(10, 1000, eigvalsSNR=100) # 1000 SPD Matrices of size 10x10 # (3)
using Plots
heatmap(Matrix(ℛ[1]), yflip=true, c=:bluesreds)
ℋ=randP(ComplexF64, 20, 1000) # 1000 Hermitian Matrices of size 20x20 # (4)
```

"""
randPosDefMat(n::Int;
              df::Int=2,
              eigvalsSNR::Real=10e3) =
    congruence(randU(n), randΛ(n, df=df, eigvalsSNR=eigvalsSNR), ℍ)


randPosDefMat(::Type{Complex{T}}, n::Int;
             df::Int=2,
             eigvalsSNR::Real=10e3) where {T<:AbstractFloat} =
    congruence(randU(ComplexF64, n), randΛ(n, df=df, eigvalsSNR=eigvalsSNR), ℍ)


function randPosDefMat(n::Int, k::Int;
             df::Int=2,
             eigvalsSNR::Real=10e3,
             SNR::Real=100,
             commuting=false)
    U=randU(n)
    𝐏=ℍVector(undef, k)
    φ=1/SNR
    if commuting
        @threads for j=1:k 𝐏[j]=ℍ(cong(U, randΛ(n, df=df, eigvalsSNR=eigvalsSNR), 𝕄)) end
    else
        𝐐=𝕄Vector(undef, k)
        @threads for j=1:k 𝐐[j]=randU(n) end
        @threads for j=1:k 𝐏[j]=ℍ(cong(U, randΛ(n, df=df, eigvalsSNR=eigvalsSNR), 𝕄)
                                +cong(𝐐[j], (φ*randΛ(n, df=df, eigvalsSNR=eigvalsSNR)), 𝕄)) end
    end
    return 𝐏
end


function randPosDefMat(::Type{Complex{T}}, n::Int, k::Int;
            df::Int=2,
            eigvalsSNR::Real=10e3,
            SNR::Real=100,
            commuting=false)            where {T<:AbstractFloat}
    U=randU(ComplexF64, n)
    𝐏=ℍVector(undef, k)
    φ=1/SNR
    if commuting
        @threads for j=1:k 𝐏[j]=ℍ(cong(U, randΛ(n, df=df, eigvalsSNR=eigvalsSNR), 𝕄)) end
    else
        𝐐=𝕄Vector(undef, k)
        @threads for j=1:k 𝐐[j]=randU(ComplexF64, n) end
        @threads for j=1:k 𝐏[j]=ℍ(cong(U, randΛ(n, df=df, eigvalsSNR=eigvalsSNR), 𝕄)
                            +cong(𝐐[j], (φ*randΛ(n, df=df, eigvalsSNR=eigvalsSNR)), 𝕄)) end
    end
    return 𝐏
end


randP=randPosDefMat


"""
    (1) regularize!(P::ℍ; <SNR=10e3>)
    (2) regularize!(𝐏::ℍVector; <SNR=10e3>)

 Add [white noise](https://bit.ly/2TN8472) to either
 - (1) a positive definite matrix ``P`` of size ``n⋅n``, or
 - (2) a 1d array ``𝐏`` of ``k`` positive definite matrices of size ``n⋅n``, of [ℍVector type](@ref).

 The added noise improves the matrix conditioning with respect to inversion.
 This is used to avoid numerical errors when decomposing these matrices
 or when evaluating some functions of their eigevalues such as the log.

 A constant value is added to all diagonal elements of (1) ``P``
 or (2) af all matrices in ``𝐏``,
 that is, on output:

 ``\\textrm{(1)}\\hspace{2pt}P\\leftarrow P+ηI``

 ``\\textrm{(2)}\\hspace{2pt}𝐏_i\\leftarrow 𝐏_i+ηI, \\hspace{2pt}\\textrm{for}\\hspace{2pt} i=1:k.``

 The amount of added noise ``η`` is determined by the `SNR`
 *<keyword argument>*, which by default is 10000. This is
 such that

 ``\\textrm{(1)}\\hspace{2pt}SNR=\\frac{\\displaystyle\\textrm{tr}(P)}{\\displaystyle\\textrm{tr}(ηI)}.``

 ``\\textrm{(2)}\\hspace{2pt}SNR=\\frac{\\displaystyle\\sum_{i=1}^{k}\\textrm{tr}(𝐏_i)}{\\displaystyle k\\hspace{1pt}\\textrm{tr}(ηI)}.``

 ``P`` in (1) must be flagged as Hermitian. See [typecasting matrices](@ref).

!!! note "Nota Bene"
    The keyword argument `SNR` expresses a SNR ([signal-to-noise ratio](https://bit.ly/1VvpvnQ)),
    and is not expressed in decibels,  but as the SNR variance ratio.
    It must be a positive number. Differently from function
    `randΛ`[`randEigvalsMat`](@ref), `randλ`[`randEigvals`](@ref) and
    `randP`[`randPosDefMat`](@ref), the SNR here is not the expected SNR,
    but the actual SNR.

**See also**: `randP` ([`randPosDefMat`](@ref)).

 **Examples**
```julia
# (1)
using LinearAlgebra, Plots, PosDefManifold
n=3
U=randU(n)
# in Q we will write two matrices,
# the unregularized and regularized matrix side by side
Q=Matrix{Float64}(undef, n, n*2)
P=ℍ(U*Diagonal(randn(n).^2)*U') # generate a real 3x3 positive matrix
for i=1:n, j=1:n Q[i, j]=P[i, j] end
regularize!(P, SNR=5)
for i=1:n, j=1:n Q[i, j+n]=P[i, j] end # the regularized matrix is on the right
heatmap(Matrix(Q), yflip=true, c=:bluesreds)

# (2)
𝐏=[ℍ(U*Diagonal(randn(3).^2)*U') for i=1:5] # 5 real 3x3 positive matrices
regularize!(𝐏, SNR=1000)

## Run a test
using LinearAlgebra
𝐏=randP(10, 100, SNR=1000); # 100 real Hermitian matrices
signalVar=sum(tr(P) for P in 𝐏);
regularize!(𝐏, SNR=1000);
signalPlusNoiseVar=sum(tr(P) for P in 𝐏);
output_snr=signalVar/(signalPlusNoiseVar-signalVar)
# output_snr should be approx. equal to 1000
```

"""
function regularize!(P::ℍ;
                    SNR=10e3)
    n=size(P, 1)
    η=tr(P)/(SNR*n)
    for i in 1:n P[i, i]+=η  end
end

function regularize!(𝐏::ℍVector;
                    SNR=10e3)
    k=length(𝐏)
    n=size(𝐏[1], 1)
    η=sum(tr(P) for P in 𝐏)/(SNR*n*k)
    for l in 1:k, i in 1:n 𝐏[l][i, i]+=η  end
end


"""
    gram(X::Matrix{T}) where T<:RealOrComplex

 Given a generic data matrix ``X``, comprised of real or complex elements,
 return the normalized [Gram matrix](https://bit.ly/2I0FQn2), that is,
 the covariance matrix of ``X``
 corrected by sample size, but without subtracting the mean.

 The result is flagged as `Hermitian`.
 See [typecasting matrices](@ref).

!!! note "Nota Bene"
    If ``X`` is wide or square (r<=c) return ``XX^H/c``.
    If ``X`` is tall (r>c)            return ``X^HX/r``.

 **Examples**
```julia
using PosDefManifold
X=randn(5, 150);
G=gram(X) # => G=X*X'/150
X=randn(100, 2);
F=gram(X); # => G=X'*X/100
```

"""
function gram(X::𝕄{T}) where T<:RealOrComplex
    (r, c)=size(X)
    r<c ? ℍ((X*X')/c) : ℍ((X'*X)/r)
end # function gram


"""
trade(P::ℍ{T}) where T<:RealOrComplex

 Given a positive definite matrix `P`, return as a 2-tuple the
 *trace* and the *determinant* of `P`.
 This is used to plot positive matrices in two dimensions
 (*TraDe plots*: log(trace/n) vs. log(determinant), see exemple here below).

 `P` must be flagged by julia as `Hermitian`.
  See [typecasting matrices](@ref).

 **Examples**
```julia
using PosDefManifold
P=randP(3)
t, d=trade(P)  # equivalent to (t, d)=trade(P)

# TraDe plot
using Plots
k=100
n=10
𝐏=randP(n, k, SNR=1000); # 100 real Hermitian matrices
x=Vector{Float64}(undef, k)
y=Vector{Float64}(undef, k)
for i=1:k
    x[i], y[i] = trade(𝐏[i])
end
x=log.(x./n)
y=log.(y)
plot(x, y, seriestype=:scatter)
```

"""
trade(P::ℍ{T}) where T<:RealOrComplex = (tr(P), det(P))
