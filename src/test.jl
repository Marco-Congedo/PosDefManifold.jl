# Unit test.jl, part of PosDefManifold Package for julia language
# v 0.1.0 - last update 1th of April 2019
#
# MIT License
# Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
# https://sites.google.com/site/marcocongedo/home
#
# DESCRIPTION
# This Unit tests all functions in PosDefManifold.
# Unce you ran it, for each method of each function,
# a star is printed if the test is succesful,
# a no_entry sign is printed if the test is not succesful.
# If there are fails, the concerned functions will be listed as Warnings
# and returned by the testall() function as an array of strings

# RUN TESTS
# just hit CRTL+A (seect all) and Shift+Enter (run),
# then invoke 'testall()' in the REPL

# USE IT FOR TESTING YOUR OWN PACKAGES
# - Copy this unit in the same directory where the Main
#     module of your library is.
# - Replace 'PosDefManifold' by the name of your module(s) in the
#    'using' command here below.
# - Overwrite the tests() function with your own tests.

push!(LOAD_PATH, pwd())
using LinearAlgebra, Statistics, PosDefManifold

function tests();
    n=10
    m=20
    P=randP(n)      # an SPD Hermitian matrix
    Q=randP(n)      # an SPD Hermitian matrix
    PC=randP(ComplexF64, n) # a complex Hermitian matrix
    T=randn(m, n)   # a real tall Matrix
    TC=randn(ComplexF64, m, n)   # a complex tall Matrix
    W=randn(n, m)   # a wide Matrix
    X=randn(n, n)   # a square matrix
    XC=randn(ComplexF64, n, n) # a complex square matrix
    λ=[0.0, 0.1, 0.2, 0.3] # a real vector
    λc=[0.0+0.0im, 0.1+0.1im, 0.2+0.2im, 0.3+0.3im] # a complex vector
    Λ=Diagonal(λ)   # a real Diagonal matrix
    ΛC=Diagonal(λc) # a complex Diagonal matrix

    P_=ℍ([  3. 2. 1.;
            2. 4. 2.;
            1. 2. 5.])
    PC_=ℍ([3.0+0.0im 2.0-0.2im 1.0-0.1im;
            2.0+0.2im 4.0+0.0im 2.0-0.2im;
            1.0+0.1im 2.0+0.2im 5.0+0.0im])

    Q_=ℍ(P_+[ 2. 1. 0.;
            1. 3. 1.;
            0. 1. 4.])
    QC_=ℍ(PC_+[2.0+0.0im 1.0-0.2im 0.0-0.1im;
            1.0+0.2im 3.0+0.0im 1.0-0.2im;
            0.0+0.1im 1.0+0.2im 4.0+0.0im])


    PC_small=randP(ComplexF64, 3)

    A3x2=[4. 3.; 2. 5.; 1. 2.]

    # functions in LinearAlgebrainP.jl
    print("Testing functions in unit 'LinearAlgebrainP.jl'")

    ## Matrix Normalizations

    name="function det1"; newTest(name)
    det(det1(P))  ≈ 1 ?  OK() : OH(name*" real case")
    det(det1(PC)) ≈ 1 ?  OK() : OH(name*" complex case")

    name="function tr1"; newTest(name)
    tr(tr1(P))  ≈ 1 ?    OK() : OH(name*" real case")
    tr(tr1(PC)) ≈ 1 ?    OK() : OH(name*" complex case")

    name="function normalizeCol!"; newTest(name)
    j=rand(1:n)
    normalizeCol!(T, j)
    norm(T[:, j])       ≈ 1 ? OK() : OH(name*" Method 1 real case")
    normalizeCol!(T, j, 2)
    norm(T[:, j])       ≈ 0.5 ? OK() : OH(name*" Method 2 real case")
    normalizeCol!(T, 2:3)
    norm(T[:, 2])+norm(T[:, 3]) ≈ 2 ? OK() : OH(name*" Method 3 real case")
    normalizeCol!(T, 2:3, 2)
    norm(T[:, 2])+norm(T[:, 3]) ≈ 1 ? OK() : OH(name*" Method 4 real case")
    normalizeCol!(TC, j)
    norm(TC[:, j])       ≈ 1 ? OK() : OH(name*" Method 1 complex case")
    normalizeCol!(TC, j, 2)
    norm(TC[:, j])       ≈ 0.5 ? OK() : OH(name*" Method 2 complex case")
    normalizeCol!(TC, 2:3)
    norm(TC[:, 2])+norm(TC[:, 3]) ≈ 2 ? OK() : OH(name*" Method 3 complex case")
    normalizeCol!(TC, 2:3, 2)
    norm(TC[:, 2])+norm(TC[:, 3]) ≈ 1 ? OK() : OH(name*" Method 4 complex case")

    ## Boolean functions of matrices

    name="function ispos"; newTest(name)
    ispos(λ, 🔔=false) == false ? OK() : OH(name*" Method 1 real case")
    ispos(Λ, 🔔=false) == false ? OK() : OH(name*" Method 2 real case")

    ## Scalar Functions of Matrices

    name="function colProd"; newTest(name)
    j1=1; j2=rand(2:n);
    s=𝚺(T[:, j1].*T[:, j2])
    colProd(T, j1, j2) ≈ s ? OK() : OH(name*" real case")
    s=𝚺(conj(TC[:, j1]).*TC[:, j2])
    colProd(TC, j1, j2) ≈ s ? OK() : OH(name*" complex case")

    name="function sumOfSqr"; newTest(name)
    sumOfSqr(P_)        ≈ 68    ? OK() : OH(name*" Method 1 real case")
    sumOfSqr(P_, 2)     ≈ 24    ? OK() : OH(name*" Method 2 real case")
    sumOfSqr(P_, 1:2)   ≈ 38    ? OK() : OH(name*" Method 3 real case")
    sumOfSqr(PC_)       ≈ 68.18 ? OK() : OH(name*" Method 1 complex case")
    sumOfSqr(PC_, 2)    ≈ 24.08 ? OK() : OH(name*" Method 2 complex case")
    sumOfSqr(PC_, 1:3)  ≈ 68.18 ? OK() : OH(name*" Method 3 complex case")

    name="function sumOfSqrDiag"; newTest(name)
    sumOfSqrDiag(P_)    ≈ 50 ? OK() : OH(name*" Method 1 real case")
    A=randn(ComplexF64, 3, 3)
    s=sum(abs2(A[i, i]) for i in 1:size(A, 1))
    sumOfSqrDiag(A)   ≈ s  ? OK() : OH(name*" Method 1 complex case")
    A=real(A)
    s=sum(A[i, i]^2 for i in 1:size(A, 1))
    sumOfSqrDiag(A)     ≈ s  ? OK() : OH(name*" Method 2")

    name="function colNorm"; newTest(name)
    colNorm(X, 2)      ≈ norm(X[:, 2]) ? OK()  : OH(name*" Method 1 real case")
    colNorm(XC, 2)     ≈ norm(XC[:, 2]) ? OK() : OH(name*" Method 1 complex case")

    name="function sumOfSqrTril"; newTest(name)
    sumOfSqrTril(A3x2, -1)≈ 9 ? OK() : OH(name*" Method 1 real case")
    A=randn(ComplexF64, 3, 2)
    s=sumOfSqr(A)-abs2(A[1, 2])
    sumOfSqrTril(A, 0)≈ s ? OK() : OH(name*" Method 1 complex case")

    name="function fidelity"; newTest(name); SKIP()
    f=fidelity(P, Q)


    ## 3. Diagonal functions of matrices

    name="function fDiagonal"; newTest(name)
    D=fDiagonal(x->x^2, P_)
    D≈Diagonal(diagm(0 => [9.,16.,25.])) ? OK() : OH(name)

    ## 4. Orthogonalization Procedures

    name="function mgs"; newTest(name)
    U=mgs(X)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=mgs(XC)
    U'*U≈I ? OK() : OH(name*" Complex Input")

    ## 5. Matrix function of matrices

    ## 6. Spectral decompositions of positive matrices

    name="function evd"; newTest(name)
    (Λ, U) = evd(P)
    U*Λ*U'≈P ? OK() : OH(name*" Real Input")
    (Λ, U) = evd(PC)
    U*Λ*U'≈PC ? OK() : OH(name*" Complex Input")


    name="function spectralFunctions"; newTest(name); SKIP()
    Q=spectralFunctions(P, x->x+1)

    name="function pow"; newTest(name)
    P½, P½ⁱ=pow(P_, 0.5, -0.5)
    P½*P½ⁱ≈I && P½*P½≈P_ ?  OK() : OH(name*" Real Input")
    P½, P½ⁱ=pow(PC_, 0.5, -0.5)
    P½*P½ⁱ≈I && P½*P½≈PC_ ? OK() : OH(name*" Complex Input")

    name="function invsqrt"; newTest(name)
    P½ⁱ=invsqrt(P_)
    P½ⁱ*P_*P½ⁱ'≈I ? OK() : OH(name)

    name="function sqr"; newTest(name)
    P²=sqr(P_)
    P² ≈ P_*P_' ? OK() : OH(name)

    name="function powerMethod"; newTest(name)
    Λ, U, iterations, covergence=powIter(ℍ(P_), size(P_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input")
    Λ, U, iterations, covergence=powIter(ℍ(PC_), size(PC_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(PC_) && U*Λ*U'≈PC_ ? OK() : OH(name*" Complex Input")

    name="function choL"; newTest(name)
    L=choL(P_)
    L*L'≈P_ ? OK() : OH(name*" Real Input")
    L=choL(PC_)
    L*L'≈PC_ ? OK() : OH(name*" Complex Input")


    # functions in SignalProcessinginP.jl
    println(" ")
    print("Testing functions in unit 'SignalProcessinginP.jl'")

    name="function randλ"; newTest(name); SKIP()
    λ=randλ(10)

    name="function randΛ"; newTest(name); SKIP()
    Λ=randΛ(10)

    name="function randU"; newTest(name);
    U=randU(10)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=randU(ComplexF64, 10)
    U'*U≈I ? OK() : OH(name*" Complex Input")

    name="function randP"; newTest(name); SKIP()

    name="function regularize!"; newTest(name)
    signalVar=tr(P)
    regularize!(P, SNR=10)
    signalPlusNoiseVar=tr(P)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 1")

    signalVar=tr(PC)
    regularize!(PC, SNR=10)
    signalPlusNoiseVar=tr(PC)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Complex Input Method 1")

    ℘=randP(5, 20)
    signalVar=𝚺(tr(P) for P in ℘)
    regularize!(℘, SNR=10)
    signalPlusNoiseVar=𝚺(tr(P) for P in ℘)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")

    ℘=randP(ComplexF64, 5, 20)
    signalVar=𝚺(tr(P) for P in ℘)
    regularize!(℘, SNR=10)
    signalPlusNoiseVar=𝚺(tr(P) for P in ℘)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")

    name="function gram"; newTest(name); SKIP()
    gramMat=gram(T)

    name="function trade"; newTest(name); SKIP()
    trade1, trade2=trade(P)

    # functions in RiemannianGeometryinP.jl
    println(" ")
    print("Testing functions in unit 'RiemannianGeometryinP.jl'")

    name="function geodesic"; newTest(name); SKIP()
    geodesic(P, Q, 0.5, Euclidean)

    name="function distanceSqr"; newTest(name); SKIP()
    d=distanceSqr(P, Wasserstein)

    name="function distance"; newTest(name); SKIP()
    d=distance(P, logdet0)

    name="function distanceSqrMatrix"; newTest(name); SKIP()

    name="function distanceMatrix"; newTest(name); SKIP()

    name="function laplacian"; newTest(name); SKIP()

    name="function laplacianEigenMaps"; newTest(name); SKIP()

    name="function spectralEmbedding"; newTest(name); SKIP()


    name="function generalizedMean"; newTest(name);
    ℘=ℍVector([P_, Q_])
    w=[0.2, 0.8]
    p=0.5
    ℍ( (P_^p+Q_^p)/2) ^(1/p) ≈ generalizedMean(℘, p) ? OK() : OH(name*" Real Input 1")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Real Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w) ? OK() : OH(name*" Real Input 3")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≉ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Real Input 4")
    ℍ( (ℍ(0.4*P_^p)+ℍ(1.6*Q_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Real Input 5")
    ℘=ℍVector([PC_, QC_])
    w=[0.2, 0.8]
    ℍ( (PC_^p+QC_^p)/2) ^(1/p) ≈ generalizedMean(℘, p) ? OK() : OH(name*" Complex Input 1")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Complex Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w) ? OK() : OH(name*" Complex Input 3")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≉ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Complex Input 4")
    ℍ( (ℍ(0.4*PC_^p)+ℍ(1.6*QC_^p))  )^(1/p) ≈ generalizedMean(℘, p, w=w, ✓w=false) ? OK() : OH(name*" Complex Input 5")


    name="function logdet0Mean"; newTest(name);
    ℘=ℍVector([P_, Q_])
    w=[0.5, 0.5]
    P½, P½ⁱ=pow(P_, 0.5, -0.5)
    GM=P½*(P½ⁱ*Q_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(℘) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 1")
    ldG, iter, conv = logdet0Mean(℘, w=w) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 2")
    ℘=℘=ℍVector([PC_, QC_])
    P½, P½ⁱ=pow(PC_, 0.5, -0.5)
    GM=P½*(P½ⁱ*QC_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(℘) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 1")
    ldG, iter, conv = logdet0Mean(℘, w=w) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 2")


end

global failed=[]

function newTest(name::String)
    sleep(0.025)
    println(" ")
    print(rpad(name*":", 30))
end

OK()=print("⭐ ")

function OH(name::String)
    print("⛔ ")
    push!(failed, name)
end

SKIP()=print(" skypped")

function testall()
    println("\n⭐ "," PosDefManifold testing utility", "⭐\n")
    println("Starting tests...")
    tests()
    # print out the tests that have failed (if any)
    println("\n")
    length(failed)==0 ? @info("All tests were succesful!") :
        for s in failed @warn("Test of $s failed") end
        #for i=1:length(failed) @warn("Test of ", failed[i]," failed") end
    if length(failed)>0
        return failed
        println("\nList of test that failed:")
    end
end # function testall

#clipboard("testall()")
#@info("\nhit CTRL+V+ENTER on the REPL for running the tests of PosDefManifold.")
#fails=testAll()
