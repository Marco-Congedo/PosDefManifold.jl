#    Main Module of the  PosDefManifold Package for julia language
#    v 0.1.2 - last update 22th of April 2019
#
#    MIT License
#    Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home


# __precompile__()

module PosDefManifold

using LinearAlgebra, Statistics

# constants
const sqrt2=√2
const invsqrt2=1/sqrt2
const minpos=1e-15
const maxpos=1e15
#const ghostMat=Matrix{Nothing}(undef, 1, 1)

# aliases
𝚺 = sum          # alias for sum, tab-completition: \bfSigma
𝛍 = mean         # alias for mean, tab-completition: \bfmu
⋱ = Diagonal     # alias for Diagonal, tab-completition: ⋱\ddots
ℍ = Hermitian    # alias for Hermitian, tab completion \bbH

# types

RealOrComplex=Union{Real, Complex}
ℍVector=Vector{ℍ}           # vector of Hermitian matrices
ℍVector₂=Vector{ℍVector}    # vector of vectors of Hermitian matrices


@enum Metric begin
    Euclidean    =1  # distance: δ_e; mean: Arithmetic
    invEuclidean =2  # distance: δ_i; mean: Harmonic
    ChoEuclidean =3  # distance: δ_c; mean: Cholesky Euclidean
    logEuclidean =4  # distance: δ_l; mean: Log Euclidean
    logCholesky  =5  # distance: δ_c; mean: Log-Cholesky
    Fisher       =6  # distance: δ_f; mean: Fisher (Cartan, Karcher, Pusz-Woronowicz,...)
    logdet0      =7  # distance: δ_s; mean: LogDet (S, α, Bhattacharyya, Jensen,...)
    Jeffrey      =8  # distance: δ_j; mean: Jeffrey (symmetrized Kullback-Leibler)
    VonNeumann   =9  # distance: δ_v; mean: Not Availale
    Wasserstein  =10 # distance: δ_w; mean: Wasserstein (Bures, Hellinger, ...)
    #...
end

import  Statistics.mean,
        LinearAlgebra.tr

export
    # From this module

    #constants
    sqrt2,
    invsqrt2,
    minpos,
    maxpos,
    ghostMat,

    #aliases
    𝚺,
    𝛍,
    ⋱,
    ℍ,

    #types
    RealOrComplex,
    ℍVector,
    ℍVector₂,
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
    det1,
    tr1,
    normalizeCol!,
    colProd,
    sumOfSqr,
    sumOfSqrDiag,
    colNorm,
    sumOfSqrTril,
    fidelity,
    fDiagonal, 𝑓𝑫,
    mgs,
    evd,
    spectralFunctions,
    ispos,
    pow,
    invsqrt,
    sqr,
    powerIterations, powIter,
    choL,

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
    distanceSqr, distance²,
    distance,
    geodesic,
    distanceSqrMat, distance²Mat,
    distanceMat,
    laplacian,
    laplacianEigenMaps, laplacianEM,
    spectralEmbedding,
    mean,
    means,
    generalizedMean,
    powerMean,
    logdet0Mean,
    wasMean,
    logMap,
    expMap,
    vecP,
    matP,
    procrustes,

    # from Test.jl
    testall


include("linearAlgebra.jl")
include("signalProcessing.jl")
include("riemannianGeometry.jl")
include("test.jl")

println("\n⭐ "," Welcome to the PosDefManifold package", " ⭐\n")

end # module
