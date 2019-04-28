#    Unit test.jl, part of PosDefManifold Package for julia language
#    v 0.1.3 - last update 28th of April 2019
#
#    MIT License
#    Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home
#
#    DESCRIPTION
#    This Unit tests all functions in PosDefManifold.
#    Unce you ran it, for each method of each function,
#    a star sign is printed if the test is succesful,
#    a play sign is printed if the function executes correctly and
#    a no_entry sign is printed if the test is not succesful.
#    If there are fails, the concerned functions will be listed as Warnings
#    and returned by the testall() function as an array of strings

#    RUN TESTS
#    just hit CRTL+A (seect all) and Shift+Enter (run),
#    then invoke 'testall()' in the REPL

#    USE IT FOR TESTING YOUR OWN PACKAGES
#    - Copy this unit in the same directory where the Main
#      module of your library is.
#    - Replace 'PosDefManifold' by the name of your module(s) in the
#      'using' command here below.
#    - Overwrite the tests() function with your own tests.

push!(LOAD_PATH, pwd())
using LinearAlgebra, Statistics, PosDefManifold

function tests();
    metrics=[Metric(i) for i in 1:10]
    n=10
    m=20
    P=randP(n)      # an SPD Hermitian matrix
    Q=randP(n)      # an SPD Hermitian matrix
    PC=randP(ComplexF64, n) # a complex Hermitian matrix
    QC=randP(ComplexF64, n) # a complex Hermitian matrix
    T=randn(m, n)   # a real tall Matrix
    TC=randn(ComplexF64, m, n)   # a complex tall Matrix
    T2=randn(m, n)   # a real tall Matrix
    TC2=randn(ComplexF64, m, n)   # a complex tall Matrix
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

    𝐏=randP(10, 4)
    𝐏C=randP(ComplexF64, 10, 4)
    𝐐=randP(10, 4)
    𝐐C=randP(ComplexF64, 10, 4)

    # functions in LinearAlgebrainP.jl
    print("Testing functions in unit 'LinearAlgebrainP.jl'")

    ## 1. Matrix Normalizations

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

    ## 2. Boolean functions of matrices

    name="function ispos"; newTest(name)
    ispos(λ, 🔔=false) == false ? OK() : OH(name*" Method 1 real case")
    ispos(Λ, 🔔=false) == false ? OK() : OH(name*" Method 2 real case")

    ## 3. Scalar Functions of Matrices

    name="function colProd"; newTest(name)
    j1=1; j2=rand(2:n);
    s=𝚺(T[:, j1].*T[:, j2])
    colProd(T, j1, j2) ≈ s ? OK() : OH(name*" Method 1 real case")
    s=𝚺(conj(TC[:, j1]).*TC[:, j2])
    colProd(TC, j1, j2) ≈ s ? OK() : OH(name*" Method 1 complex case")
    s=𝚺(T[:, j1].*T2[:, j2])
    colProd(T, T2, j1, j2) ≈ s ? OK() : OH(name*" Method 2 real case")
    s=𝚺(conj(TC[:, j1]).*TC2[:, j2])
    colProd(TC, TC2, j1, j2) ≈ s ? OK() : OH(name*" Method 2 complex case")


    name="function sumOfSqr"; newTest(name)
    sumOfSqr(Matrix(P_))  ≈ 68    ? OK() : OH(name*" Method 1 real case")
    sumOfSqr(P_)          ≈ 68    ? OK() : OH(name*" Method 2 real case")
    sumOfSqr(P_, 2)       ≈ 24    ? OK() : OH(name*" Method 3 real case")
    sumOfSqr(P_, 1:2)     ≈ 38    ? OK() : OH(name*" Method 4 real case")
    sumOfSqr(Matrix(PC_)) ≈ 68.18 ? OK() : OH(name*" Method 1 complex case")
    sumOfSqr(PC_)         ≈ 68.18 ? OK() : OH(name*" Method 2 complex case")
    sumOfSqr(PC_, 2)      ≈ 24.08 ? OK() : OH(name*" Method 3 complex case")
    sumOfSqr(PC_, 1:3)    ≈ 68.18 ? OK() : OH(name*" Method 4 complex case")


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


    name="function tr"; newTest(name)
    tr(P, Q) ≈ tr(P*Q) ? OK() : OH(name*" Method 1 real case")
    tr(PC, QC) ≈ tr(PC*QC) ? OK() : OH(name*" Method 1 complex case")
    tr(P, X) ≈ tr(P*X) ? OK() : OH(name*" Method 2 real case")
    tr(PC, XC) ≈ tr(PC*XC) ? OK() : OH(name*" Method 2 complex case")


    name="function quadraticForm"; newTest(name)
    v=randn(n)
    vC=randn(ComplexF64, n)
    quadraticForm(v, Matrix(P)) ≈ v'*P*v ? OK() : OH(name*" Method 1 real case")
    quadraticForm(v, P) ≈ v'*P*v ? OK() : OH(name*" Method 2 real case")
    quadraticForm(v, LowerTriangular(Matrix(P))) ≈ v'*P*v ? OK() : OH(name*" Method 3 real case")
    quadraticForm(vC, Matrix(PC)) ≈ vC'*PC*vC ? OK() : OH(name*" Method 1 complex case")
    quadraticForm(vC, PC) ≈ vC'*PC*vC ? OK() : OH(name*" Method 2 complex case")



    name="function fidelity"; newTest(name);
    # Test compilation only
    f=fidelity(P, Q); RUN()
    f=fidelity(PC, QC); RUN()


    ## 4. Diagonal functions of matrices

    name="function fDiagonal"; newTest(name)
    D=fDiagonal(x->x^2, P_)
    D≈Diagonal(diagm(0 => [9.,16.,25.])) ? OK() : OH(name)

    ## 5. Unitary functions of matrices

    name="function mgs"; newTest(name)
    U=mgs(X)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=mgs(XC)
    U'*U≈I ? OK() : OH(name*" Complex Input")

    ## 6. Matrix function of matrices

    ## 7. Spectral decompositions of positive matrices

    name="function evd"; newTest(name)
    (Λ, U) = evd(P)
    U*Λ*U'≈P ? OK() : OH(name*" Real Input")
    (Λ, U) = evd(PC)
    U*Λ*U'≈PC ? OK() : OH(name*" Complex Input")


    name="function spectralFunctions"; newTest(name);
    spectralFunctions(P, x->x+1); RUN()
    spectralFunctions(PC, abs2); RUN()


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


    name="function powerIterations"; newTest(name)
    Λ, U, iterations, convergence=powIter(P_, size(P_, 2), evalues=true; tol=tol=1e-9)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 1")
    Λ, U, iterations, convergence=powIter(PC_, size(PC_, 2), evalues=true; tol=tol=1e-9)
    sort(diag(Λ))≈eigvals(PC_) && U*Λ*U'≈PC_ ? OK() : OH(name*" Complex Input Method 1")
    Λ, U, iterations, convergence=powIter(ℍ(P_), size(P_, 2), evalues=true; tol=tol=1e-9)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 2")
    Λ, U, iterations, convergence=powIter(ℍ(PC_), size(PC_, 2), evalues=true; tol=tol=1e-9)
    sort(diag(Λ))≈eigvals(PC_) && U*Λ*U'≈PC_ ? OK() : OH(name*" Complex Input Method 2")
    Λ, U, iterations, convergence=powIter(LowerTriangular(P_), size(P_, 2), evalues=true; tol=tol=1e-9)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 3")


    name="function choL"; newTest(name)
    L=choL(P_)
    L*L'≈P_ ? OK() : OH(name*" Real Input")
    L=choL(PC_)
    L*L'≈PC_ ? OK() : OH(name*" Complex Input")

    # 8. Decompositions involving triangular matrices

    # functions in SignalProcessinginP.jl
    println(" ")
    print("Testing functions in unit 'SignalProcessinginP.jl'")

    name="function randλ"; newTest(name);
    randλ(10); RUN()


    name="function randΛ"; newTest(name);
    randΛ(10); RUN()


    name="function randU"; newTest(name);
    U=randU(10)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=randU(ComplexF64, 10)
    U'*U≈I ? OK() : OH(name*" Complex Input")


    name="function randP"; newTest(name);
    randP(3); RUN()
    randP(ComplexF64, 3); RUN()


    name="function regularize!"; newTest(name)
    signalVar=tr(P)
    regularize!(P, SNR=10)
    signalPlusNoiseVar=tr(P)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 1")

    signalVar=tr(PC)
    regularize!(PC, SNR=10)
    signalPlusNoiseVar=tr(PC)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Complex Input Method 1")

    𝐏=randP(5, 20)
    signalVar=𝚺(tr(P) for P in 𝐏)
    regularize!(𝐏, SNR=10)
    signalPlusNoiseVar=𝚺(tr(P) for P in 𝐏)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")

    𝐏=randP(ComplexF64, 5, 20)
    signalVar=𝚺(tr(P) for P in 𝐏)
    regularize!(𝐏, SNR=10)
    signalPlusNoiseVar=𝚺(tr(P) for P in 𝐏)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")


    name="function gram"; newTest(name);
    gram(T); RUN()
    gram(W); RUN()


    name="function trade"; newTest(name);
    trade(P); RUN()
    trade(PC); RUN()


    # functions in RiemannianGeometryinP.jl
    println(" ")
    print("Testing functions in unit 'RiemannianGeometryinP.jl'")

    name="function geodesic"; newTest(name);
    (geodesic(m, P, Q, 0.5) for m in metrics if m≠9); RUN()
    (geodesic(m, PC, QC, 0.5) for m in metrics if m≠9); RUN()


    name="function distanceSqr"; newTest(name);
    for m in metrics distanceSqr(m, P) end; RUN()
    for m in metrics distanceSqr(m, P, Q) end; RUN()
    for m in metrics distanceSqr(m, PC) end; RUN()
    for m in metrics distanceSqr(m, PC, QC) end; RUN()


    name="function distance"; newTest(name);
    # since it calls distanceSqr just check the call works
    for m in metrics distance(m, P) end; RUN()
    for m in metrics distance(m, P, Q) end; RUN()
    for m in metrics distance(m, PC) end; RUN()
    for m in metrics d=distance(m, PC, QC); end; RUN()


    name="function distanceSqrMat (I)"; newTest(name);
    k=length(𝐏)
    for m in metrics
            D=distanceSqrMat(m, 𝐏)
            manualD=LowerTriangular(Matrix{Float64}(undef, k, k))
            for j=1:k, i=j:k manualD[i, j]=distanceSqr(m, 𝐏[i], 𝐏[j]) end
            manualD≈D ? OK() : OH(name*" Real Input, metric "*string(m))
    end


    name="function distanceSqrMat (II)"; newTest(name);
    k=length(𝐏C)
    for m in metrics
            D=distanceSqrMat(m, 𝐏C)
            manualD=LowerTriangular(Matrix{Float64}(undef, k, k))
            for j=1:k, i=j:k manualD[i, j]=distanceSqr(m, 𝐏C[i], 𝐏C[j]) end
            manualD≈D ? OK() : OH(name*" Complex Input, metric "*string(m))
    end


    name="function distanceMat"; newTest(name); SKIP()


    name="function laplacian"; newTest(name);
    Dsqr=distanceSqrMat(logEuclidean, 𝐏)
    lap=laplacian(Dsqr); RUN()


    name="function laplacianEigenMaps"; newTest(name);
    laplacianEM(lap, 2); RUN()


    name="function spectralEmbedding"; newTest(name);
    spectralEmbedding(logEuclidean, 𝐏, 2); RUN()


    name="function mean"; newTest(name);
    (mean(m, P, Q) for m in metrics); RUN()
    (mean(m, 𝐏) for m in metrics); RUN()
    (mean(m, PC, QC) for m in metrics); RUN()
    (mean(m, 𝐏C) for m in metrics); RUN()


    name="function means"; newTest(name);
    means(logEuclidean, ℍVector₂([𝐏, 𝐐])); RUN()
    means(logEuclidean, ℍVector₂([𝐏C, 𝐐C])); RUN()


    name="function generalizedMean"; newTest(name);
    𝐏=ℍVector([P_, Q_])
    w=[0.2, 0.8]
    p=0.5
    ℍ( (P_^p+Q_^p)/2) ^(1/p) ≈ generalizedMean(𝐏, p) ? OK() : OH(name*" Real Input 1")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w) ? OK() : OH(name*" Real Input 3")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≉ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 4")
    ℍ( (ℍ(0.4*P_^p)+ℍ(1.6*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 5")
    𝐏=ℍVector([PC_, QC_])
    w=[0.2, 0.8]
    ℍ( (PC_^p+QC_^p)/2) ^(1/p) ≈ generalizedMean(𝐏, p) ? OK() : OH(name*" Complex Input 1")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w) ? OK() : OH(name*" Complex Input 3")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≉ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 4")
    ℍ( (ℍ(0.4*PC_^p)+ℍ(1.6*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 5")


    name="function logdet0Mean"; newTest(name);
    w=[0.5, 0.5]
    P½, P½ⁱ=pow(P_, 0.5, -0.5)
    GM=P½*(P½ⁱ*Q_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(ℍVector([P_, Q_])) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 1")
    ldG, iter, conv = logdet0Mean(ℍVector([P_, Q_]); w=w) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 2")
    P½, P½ⁱ=pow(PC_, 0.5, -0.5)
    GM=P½*(P½ⁱ*QC_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(ℍVector([PC_, QC_])) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 1")
    ldG, iter, conv = logdet0Mean(ℍVector([PC_, QC_]); w=w) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 2")


    name="function wasMean"; newTest(name);
    wasMean(ℍVector([P_, Q_])); RUN()
    wasMean(ℍVector([PC_, QC_])); RUN()


    name="function powerMean"; newTest(name);
    powerMean(ℍVector([P_, Q_]), 0.5); RUN()
    powerMean(ℍVector([PC_, QC_]), 0.5); RUN()


    name="function logMap"; newTest(name);
    logMap(Fisher, P, Q); RUN()
    logMap(Fisher, PC, QC); RUN()


    name="function expMap"; newTest(name);
    expMap(Fisher, P, Q); RUN()
    expMap(Fisher, PC, QC); RUN()

    G=mean(Fisher, P, Q)
    S=logMap(Fisher, P, G)
    H=expMap(Fisher, S, G)
    if P≉H
        println("")
        @error("either logMap or expMap or both do not give the expected output in the real case")
    end

    G=mean(Fisher, PC, QC)
    S=logMap(Fisher, PC, G)
    H=expMap(Fisher, S, G)
    if PC≉H
        println("")
        @error("either logMap or expMap or both do not give the expected output in the complex case")
    end


    name="function vecP"; newTest(name);
    v=vecP(P); RUN()
    vC=vecP(PC); RUN()


    name="function matP"; newTest(name);
    Pnew=matP(v); RUN()
    PCnew=matP(vC); RUN()

    if P≉Pnew
        println("")
        @warn("either vecP or matP or both do not give the expected output in the real case")
    end

    if PC≉PCnew
        println("")
        @warn("either vecP or matP or both do not give the expected output in the complex case")
    end

    name="function procrustes"; newTest(name);
    procrustes(P, Q); RUN()

end # function tests

function newTest(name::String)
    sleep(0.025)
    println(" ")
    print(rpad(name*":", 30))
end

OK()=print("⭐ ")

RUN()=print("▶ ")

function OH(name::String)
    print("⛔ ")
    push!(failing_tests, name)
end

SKIP()=print(" skypped")

failing_tests=[]

function testall()
    println("\n⭐ "," PosDefManifold testing utility", "⭐\n")
    println("Starting tests...")
    tests()
    # print out the tests that have failed (if any)
    println("\n")
    length(failing_tests)==0 ? @info("All tests were succesful!") :
        for s in failing_tests @warn("Test of $s failed") end
    for i=1:length(failing_tests) pop!(failing_tests) end
    return nothing
end # function testall

#clipboard("testall()")
#@info("\nhit CTRL+V+ENTER on the REPL for running the tests of PosDefManifold.")
#fails=testAll()
