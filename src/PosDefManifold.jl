#   Main Module of the  PosDefManifold Package for julia language
#   v0.5.2 - May 2025

#   MIT License
#   Copyright (c) 2019-25, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home

# __precompile__()

module PosDefManifold

using LinearAlgebra, Statistics, Base.Threads

# Special instructions and variables
BLAS.set_num_threads(Sys.CPU_THREADS)

# constants
const sqrt2     =   √2
const invsqrt2  =   1/sqrt2
const golden    =   (√5+1)/2 # 1.618033988749894848204586834365638117720309...
const goldeninv =   golden-1.0 # = 1/gold = (√5-1)/2 = ...
const maxpos    =   1e15

const titleFont     = "\x1b[91m " # light red
const separatorFont = "\x1b[31m " # red
const defaultFont   = "\x1b[0m"
const greyFont      = "\x1b[90m"
const 📌            = titleFont*"PosDefManifold.jl "*defaultFont

# aliases
𝕄 = Matrix	            # tab-completion: \bbM
𝔻 = Diagonal	        # tab-completition: \bbD
ℍ = Hermitian          # tab completion \bbH
𝕃 = LowerTriangular    # tab completition \bbL

# private aliases
𝚺 = sum                 # tab-completition: \bfSigma
𝛍 = mean                # tab-completition: \bfmu

# types
MatrixVector = 𝕄Vector = Vector{𝕄}            # vector of Matrices
MatrixVector₂= 𝕄Vector₂= Vector{𝕄Vector}      # vector of vectors of Matrices

DiagonalVector = 𝔻Vector = Vector{𝔻}          # vector of Diagonal Matrices
DiagonalVector₂= 𝔻Vector₂= Vector{𝔻Vector}    # vector of vectors of Diagonal matrices

LowerTriangularVector = 𝕃Vector = Vector{𝕃}       # vector of LowerTriangular Matrices
LowerTriangularVector₂= 𝕃Vector₂= Vector{𝕃Vector} # vector of vectors of Diagonal matrices

HermitianVector = ℍVector = Vector{ℍ}           # vector of Hermitian matrices
HermitianVector₂= ℍVector₂= Vector{ℍVector}     # vector of vectors of Hermitian matrices

RealOrComplex=Union{Real, Complex}

#AnyMatrix=Union{𝔻{T}, 𝕃{T}, ℍ{T}, 𝕄{T}} where T<:RealOrComplex
AnyMatrix=Union{𝔻, 𝕃, ℍ, 𝕄}
AnyMatrixVector=Union{𝕄Vector, 𝔻Vector, 𝕃Vector, ℍVector}
AnyMatrixVector₂=Union{𝕄Vector₂, 𝔻Vector₂, 𝕃Vector₂, ℍVector₂}

# Metrics for the Manifold
@enum Metric begin
    Euclidean    =  1  # distance: δ_e; mean: Arithmetic
    invEuclidean =  2  # distance: δ_i; mean: Harmonic
    ChoEuclidean =  3  # distance: δ_c; mean: Cholesky Euclidean
    logEuclidean =  4  # distance: δ_l; mean: Log Euclidean
    logCholesky  =  5  # distance: δ_c; mean: Log-Cholesky
    Fisher       =  6  # distance: δ_f; mean: Fisher (Cartan, Karcher, Pusz-Woronowicz,...)
    logdet0      =  7  # distance: δ_s; mean: LogDet (S, α, Bhattacharyya, Jensen,...)
    Jeffrey      =  8  # distance: δ_j; mean: Jeffrey (symmetrized Kullback-Leibler)
    VonNeumann   =  9  # distance: δ_v; mean: Not Availale
    Wasserstein  = 10 # distance: δ_w; mean: Wasserstein (Bures, Hellinger, ...)
    #...
end


import
    Statistics.mean,
    Statistics.std,
    LinearAlgebra.tr

export
    # From this module

    #constants
    sqrt2,
    invsqrt2,
    golden,
    goldeninv,
    maxpos,

    #aliases
    𝕄,
    𝔻,
    ℍ,
    𝕃,

    #types
    MatrixVector, 𝕄Vector,
    MatrixVector₂, 𝕄Vector₂,
    DiagonalVector, 𝔻Vector,
    DiagonalVector₂, 𝔻Vector₂,
    LowerTriangularVector, 𝕃Vector,
    LowerTriangularVector₂, 𝕃Vector₂,
    HermitianVector, ℍVector,
    HermitianVector₂, ℍVector₂,
    AnyMatrix,
    AnyMatrixVector,
    AnyMatrixVector₂,
    Metric,
        Euclidean,
        invEuclidean,
        ChoEuclidean,
        logEuclidean,
        logCholesky,
        Fisher,
        logdet0,
        Jeffrey,
        VonNeumann,
        Wasserstein,

    # from LinearAlgebra.jl
    typeofMatrix, typeofMat,
    typeofVector, typeofVec,
    dim,
    remove,
    isSquare,
    det1,
    tr1,
    nearestPosDef,
    nearestOrthogonal, nearestOrth,
    normalizeCol!,
    colProd,
    sumOfSqr, ss,
    sumOfSqrDiff, ssdiff,
    sumOfSqrDiag, ssd,
    colNorm,
    sumOfSqrTril, sst,
    tr,
    quadraticForm, qf,
    fidelity,
    fDiag, 𝑓𝔻,
    DiagOfProd, dop,
    mgs,
    fVec,
    congruence, cong,
    evd,
    frf,
    invfrf,
    spectralFunctions,
    ispos,
    pow,
    invsqrt,
    sqr,
    powerIterations, powIter,
    choL,
    choInv,
    choInv!,

    # from SignalProcessing.jl
    randChi², randχ²,
    randEigvals, randλ,
    randEigvalsMat, randΛ,
    randUnitaryMat, randOrthMat, randU,
    randPosDefMat, randP,
    regularize!,
    gram,
    trade,

    # from RiemannianGeometry.jl
    geodesic,
    distanceSqr, distance²,
    distance,
    distanceSqrMat, distance²Mat,
    distanceMat,
    laplacian,
    laplacianEigenMaps, laplacianEM,
    spectralEmbedding, spEmb,
    mean,
    means,
    generalizedMean,
    geometricMean, gMean,
    geometricpMean, gpMean,
    logdet0Mean, ld0Mean,
    wasMean,
    powerMean,
    inductiveMean,
    midrange,
    logMap,
    expMap,
    vecP,
    matP,
    parallelTransport, pt,
    procrustes,

    # from statistics.jl
    # mean # already exported by RiemannianGeometry.jl
    softmax,
    std,

    # from Test.jl
    testall


include("linearAlgebra.jl")
include("signalProcessing.jl")
include("riemannianGeometry.jl")
include("statistics.jl")
include("test.jl")

println("\n⭐ "," Welcome to the ", 📌, " package", " ⭐\n")
@info " "
println(" Your Machine `", separatorFont, gethostname(), defaultFont, "` (",Sys.MACHINE, ")")
println(" runs on kernel ",Sys.KERNEL," with word size ",Sys.WORD_SIZE,".")
println(" CPU  Threads: ", separatorFont, Sys.CPU_THREADS, defaultFont)
println(" Base.Threads: ", separatorFont, "$(Threads.nthreads())", defaultFont)
println(" BLAS Threads: ", separatorFont, "$(BLAS.get_num_threads())", "\n", defaultFont)

end # module
