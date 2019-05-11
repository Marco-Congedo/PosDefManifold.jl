#    Main Module of the  PosDefManifold Package for julia language
#    v 0.1.3 - last update 28th of April 2019
#
#    MIT License
#    Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home


# __precompile__()

module PosDefManifold

using LinearAlgebra, Statistics, Base.Threads

# Special instructions and variables
BLAS.set_num_threads(Sys.CPU_THREADS-Threads.nthreads())

# constants
const sqrt2=√2
const invsqrt2=1/sqrt2
const maxpos=1e15


# aliases
𝚺 = sum                 # tab-completition: \bfSigma
𝛍 = mean                # tab-completition: \bfmu
𝕄 = Matrix	            # tab-completion: \bbM
𝔻 = Diagonal	        # tab-completition: \bbD
ℍ = Hermitian          # tab completion \bbH
𝕃 = LowerTriangular    # tab completition \bbL

# types

RealOrComplex=Union{Real, Complex}

MatrixVector=Vector{𝕄}                 # vector of Matrices
𝕄Vector=MatrixVector                   # alias
MatrixVector₂=Vector{𝕄Vector}          # vector of vectors of Matrices
𝕄Vector₂=MatrixVector₂                 # alias

DiagonalVector=Vector{𝔻}               # vector of Diagonal Matrices
𝔻Vector=DiagonalVector                 # alias
DiagonalVector₂=Vector{𝔻Vector}        # vector of vectors of Diagonal matrices
𝔻Vector₂=DiagonalVector₂               # alias

LowerTriangularVector=Vector{𝕃}        # vector of LowerTriangular Matrices
𝕃Vector=LowerTriangularVector          # alias
LowerTriangularVector₂=Vector{𝕃Vector} # vector of vectors of Diagonal matrices
𝕃Vector₂=LowerTriangularVector₂        # alias

HermitianVector=Vector{ℍ}               # vector of Hermitian matrices
ℍVector=HermitianVector                 # alias
HermitianVector₂=Vector{ℍVector}        # vector of vectors of Hermitian matrices
ℍVector₂=HermitianVector₂               # alias

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


import  Statistics.mean,
        LinearAlgebra.tr

export
    # From this module

    #constants
    sqrt2,
    invsqrt2,
    minpos,
    maxpos,

    #aliases
    𝚺,
    𝛍,
    𝕄,
    𝔻,
    ℍ,
    𝕃,

    #types
    RealOrComplex,
    MatrixVector, 𝕄Vector,
    MatrixVector₂, 𝕄Vector₂,
    DiagonalVector, 𝔻Vector,
    DiagonalVector₂, 𝔻Vector₂,
    LowerTriangularVector, 𝕃Vector,
    LowerTriangularVector₂, 𝕃Vector₂,
    HermitianVector, ℍVector,
    HermitianVector₂, ℍVector₂,
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
    sumOfSqr, ss,
    sumOfSqrDiag, ssd,
    colNorm,
    sumOfSqrTril, sst,
    tr,
    quadraticForm, qf,
    fidelity,
    fDiag, 𝑓𝔻,
    DiagOfProd, dop,
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
    geodesic,
    distanceSqr, distance²,
    distance,
    distanceSqrMat, distance²Mat,
    distanceSqrMat⏩, distance²Mat⏩,
    distanceMat,
    laplacian,
    laplacianEigenMaps, laplacianEM,
    spectralEmbedding,
    mean,
    means,
    generalizedMean,
    geometricMean, gMean,
    logdet0Mean, ld0Mean,
    wasMean,
    powerMean,
    logMap,
    expMap,
    vecP,
    matP,
    procrustes,

    # from classification.jl
    softmax,

    # from Test.jl
    testall


include("linearAlgebra.jl")
include("signalProcessing.jl")
include("riemannianGeometry.jl")
include("classification.jl")
include("test.jl")

println("\n⭐ "," Welcome to the PosDefManifold package", " ⭐\n")
@info(" ")
println(" Your Machine ",gethostname()," (",Sys.MACHINE, ")")
println(" runs on kernel ",Sys.KERNEL," with word size ",Sys.WORD_SIZE,".")
println(" CPU  Threads: ",Sys.CPU_THREADS)
# Sys.BINDIR # julia bin directory
println(" Base.Threads: ", "$(Threads.nthreads())")
println(" BLAS Threads: ", "$(Sys.CPU_THREADS-Threads.nthreads())", "\n")



end # module
