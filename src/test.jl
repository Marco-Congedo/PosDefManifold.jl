#    Unit test.jl, part of PosDefManifold Package for julia language
#
#    MIT License
#    Copyright (c) 2019-25, Marco Congedo, CNRS, Grenobe, France:
#    https://sites.google.com/site/marcocongedo/home
#
#    DESCRIPTION
#    This Unit tests all functions in PosDefManifold.
#    Some functions are fully tested, the others are just executed.
#    Unce you run it, for each method of each function,
#    a ⭐ sign is printed if the test is succesful, while
#    a ⛔ sign is printed if the test is not succesful.
#    A ☆ sign is printed if the function has been executed correctly.
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

    D_=Diagonal(P_)
    DC_=Diagonal(PC_)
    E_=Diagonal(Q_)
    EC_=Diagonal(QC_)

    PC_small=randP(ComplexF64, 3)

    A3x2=[4. 3.; 2. 5.; 1. 2.]

    𝐏=randP(8, 4)
    𝐏C=randP(ComplexF64, 8, 4)
    𝐐=randP(8, 4)
    𝐐C=randP(ComplexF64, 8, 4)
    𝐃=randΛ(8, 4)
    𝐄=randΛ(8, 4)
    weights=[0.1, 0.2, 0.3, 0.4]

    # functions in LinearAlgebrainP.jl
    print("- Unit 'linearAlgebra.jl'")

    ## 1. Utilities
    name="typeofMatrix"; newTest(name)
    typeofMatrix(𝐏); RUN()

    name="dim"; newTest(name)
    dim(𝐏); RUN()

    ## 2. Matrix Normalizations and approximations

    name="det1"; newTest(name)
    det(det1(P))  ≈ 1 ?  OK() : OH(name*" real case")
    det(det1(PC)) ≈ 1 ?  OK() : OH(name*" complex case")

    name="tr1"; newTest(name)
    tr(tr1(P))  ≈ 1 ?    OK() : OH(name*" real case")
    tr(tr1(PC)) ≈ 1 ?    OK() : OH(name*" complex case")

    name="nearestPosDef"; newTest(name)
    nearestPosDef(X); RUN()
    nearestPosDef(XC); RUN()

    name="normalizeCol!"; newTest(name)
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

    ## 3. Boolean functions of matrices

    name="ispos"; newTest(name)
    ispos(λ, 🔔=false) == false ? OK() : OH(name*" Method 1 real case")
    ispos(Λ, 🔔=false) == false ? OK() : OH(name*" Method 2 real case")

    ## 4. Scalar Functions of Matrices

    name="colProd"; newTest(name)
    j1=1; j2=rand(2:n);
    s=sum(T[:, j1].*T[:, j2])
    colProd(T, j1, j2) ≈ s ? OK() : OH(name*" Method 1 real case")
    s=sum(conj(TC[:, j1]).*TC[:, j2])
    colProd(TC, j1, j2) ≈ s ? OK() : OH(name*" Method 1 complex case")
    s=sum(T[:, j1].*T2[:, j2])
    colProd(T, T2, j1, j2) ≈ s ? OK() : OH(name*" Method 2 real case")
    s=sum(conj(TC[:, j1]).*TC2[:, j2])
    colProd(TC, TC2, j1, j2) ≈ s ? OK() : OH(name*" Method 2 complex case")


    name="colNorm"; newTest(name)
    colNorm(X, 2)      ≈ norm(X[:, 2]) ? OK()  : OH(name*" Method 1 real case")
    colNorm(XC, 2)     ≈ norm(XC[:, 2]) ? OK() : OH(name*" Method 1 complex case")

    name="sumOfSqr"; newTest(name)
    sumOfSqr(Matrix(P_))  ≈ 68    ? OK() : OH(name*" Method 1 real case")
    sumOfSqr(P_)          ≈ 68    ? OK() : OH(name*" Method 2 real case")
    sumOfSqr(LowerTriangular(P_)) ≈ 59 ? OK() : OH(name*" Method 3 real case")
    sumOfSqr(Diagonal(P_)) ≈ 50 ? OK() : OH(name*" Method 4 real case")
    sumOfSqr(P_, 2)       ≈ 24    ? OK() : OH(name*" Method 5 real case")
    sumOfSqr(P_, 1:2)     ≈ 38    ? OK() : OH(name*" Method 6 real case")
    sumOfSqr(Matrix(PC_)) ≈ 68.18 ? OK() : OH(name*" Method 1 complex case")
    sumOfSqr(PC_)         ≈ 68.18 ? OK() : OH(name*" Method 2 complex case")
    sumOfSqr(LowerTriangular(PC_)) ≈ 59.089999999999996 ? OK() : OH(name*" Method 3 complex case")
    sumOfSqr(Diagonal(PC_)) ≈ 50 ? OK() : OH(name*" Method 4 complex case")
    sumOfSqr(PC_, 2)      ≈ 24.08 ? OK() : OH(name*" Method 5 complex case")
    sumOfSqr(PC_, 1:3)    ≈ 68.18 ? OK() : OH(name*" Method 6 complex case")

    name="sumOfSqrDiff"; newTest(name)
    sumOfSqrDiff(P, Q)  ≈ sumOfSqr(P-Q) ? OK() : OH(name*" Method 1 real case")
    sumOfSqrDiff(𝕃(P), 𝕃(Q))  ≈ sumOfSqr(𝕃(P)-𝕃(Q)) ? OK() : OH(name*" Method 2 real case")
    sumOfSqrDiff(𝔻(P), 𝔻(Q))  ≈ sumOfSqr(𝔻(P)-𝔻(Q)) ? OK() : OH(name*" Method 3 real case")
    sumOfSqrDiff(PC, QC)  ≈ sumOfSqr(PC-QC) ? OK() : OH(name*" Method 1 complex case")
    sumOfSqrDiff(𝕃(PC), 𝕃(QC))  ≈ sumOfSqr(𝕃(PC)-𝕃(QC)) ? OK() : OH(name*" Method 2 complex case")
    sumOfSqrDiff(𝔻(PC), 𝔻(QC))  ≈ sumOfSqr(𝔻(PC)-𝔻(QC)) ? OK() : OH(name*" Method 3 complex case")    
    
    name="sumOfSqrDiag"; newTest(name)
    sumOfSqrDiag(P_)    ≈ 50 ? OK() : OH(name*" Method 1 real case")
    A=randn(ComplexF64, 3, 3)
    s=sum(abs2(A[i, i]) for i in 1:size(A, 1))
    sumOfSqrDiag(A)   ≈ s  ? OK() : OH(name*" Method 1 complex case")
    A=real(A)
    s=sum(A[i, i]^2 for i in 1:size(A, 1))
    sumOfSqrDiag(A)     ≈ s  ? OK() : OH(name*" Method 2")


    name="sumOfSqrTril"; newTest(name)
    sumOfSqrTril(A3x2, -1)≈ 9 ? OK() : OH(name*" Method 1 real case")
    A=randn(ComplexF64, 3, 2)
    s=sumOfSqr(A)-abs2(A[1, 2])
    sumOfSqrTril(A, 0)≈ s ? OK() : OH(name*" Method 1 complex case")


    name="tr"; newTest(name)
    tr(P, Q) ≈ tr(P*Q) ? OK() : OH(name*" Method 1 real case")
    tr(PC, QC) ≈ tr(PC*QC) ? OK() : OH(name*" Method 1 complex case")
    tr(P, X) ≈ tr(P*X) ? OK() : OH(name*" Method 2 real case")
    tr(PC, XC) ≈ tr(PC*XC) ? OK() : OH(name*" Method 2 complex case")
    tr(D_, P_) ≈ tr(D_*P_) ? OK() : OH(name*" Method 3 real case")
    tr(DC_, PC_) ≈ tr(DC_*PC_) ? OK() : OH(name*" Method 3 complex case")
    tr(P_, D_) ≈ tr(D_*P_) ? OK() : OH(name*" Method 4 real case")
    tr(PC_, DC_) ≈ tr(DC_*PC_) ? OK() : OH(name*" Method 4 complex case")



    name="quadraticForm"; newTest(name)
    v=randn(n)
    vC=randn(ComplexF64, n)
    quadraticForm(v, P) ≈ v'*P*v ? OK() : OH(name*" Method 1 real case")
    quadraticForm(v, LowerTriangular(Matrix(P))) ≈ v'*Symmetric((LowerTriangular(P))')*v ? OK() : OH(name*" Method 2 real case")
    quadraticForm(v, Matrix(P), true) ≈ v'*P*v ? OK() : OH(name*" Method 3a real case")
    quadraticForm(v, Matrix(P), false) ≈ v'*P*v ? OK() : OH(name*" Method 3b real case")
    quadraticForm(vC, PC) ≈ vC'*PC*vC ? OK() : OH(name*" Method 4a complex case")
    quadraticForm(vC, Matrix(PC)) ≈ vC'*PC*vC ? OK() : OH(name*" Method 4b complex case")
    quadraticForm(vC, LowerTriangular(Matrix(PC))) ≈ vC'*LowerTriangular(Matrix(PC))*vC ? OK() : OH(name*" Method 4c complex case")

    name="fidelity"; newTest(name);
    # Test compilation only
    f=fidelity(P, Q); RUN()
    f=fidelity(PC, QC); RUN()


    ## 5. Diagonal functions of matrices

    name="fDiag"; newTest(name)
    D=fDiag(x->x^2, P_)
    D≈Diagonal(diagm(0 => [9.,16.,25.])) ? OK() : OH(name)


    name="DiagOfProd"; newTest(name)
    DiagOfProd(P, Q)≈Diagonal(P*Q) ? OK() : OH(name*" Real Input")
    DiagOfProd(PC, QC)≈Diagonal(PC*QC) ? OK() : OH(name*" Complex Input")


    ## 6. Unitary functions of matrices

    name="mgs"; newTest(name)
    U=mgs(X)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=mgs(XC)
    U'*U≈I ? OK() : OH(name*" Complex Input")

    ## 7. Matrix function of matrices

    name="fVec"; newTest(name)
    Pset=randP(4, 64)
    mean(Pset)≈fVec(mean, Pset) ? OK() : OH(name*" simple mean, real")
    mean(sqrt, Pset)≈fVec(mean, sqrt, Pset) ? OK() : OH(name*" mean, sqrt, real")
    fVecw=(randn(64)).^2
    fVecw=fVecw./sum(fVecw) # generate normalized random weights
    mean(Pset[i]*fVecw[i] for i=1:64)≈fVec(mean, Pset; w=fVecw) ? OK() : OH(name*" simple weighted mean, real")
    mean(sqrt(Pset[i])*fVecw[i] for i=1:64)≈fVec(mean, sqrt, Pset; w=fVecw) ? OK() : OH(name*" weighted mean, sqrt, real")

    name="congruence (I)"; newTest(name)
    M=randn(n, n)
    congruence(M, P, ℍ)≈M*P*M' ? OK() : OH(name)

    name="congruence (II)"; newTest(name)
    Pset=randP(n, n*4)
    Qset=cong(M, Pset, ℍVector)
    Qset≈[M*Pset[i]*M' for i=1:n*4] ? OK() : OH(name)

    name="congruence (III)"; newTest(name)
    Pset1=randP(4, 10); # generate 10 positive definite 4x4 matrices
    Pset2=randP(4, 8);
    Pset=ℍVector₂([Pset1, Pset2]);
    M=randn(4, 4)
    Qset=cong(M, Pset, MatrixVector₂)
    Qset[1][1]≈M*Pset[1][1]*M' &&
    Qset[1][5]≈M*Pset[1][5]*M' &&
    Qset[2][1]≈M*Pset[2][1]*M' &&
    Qset[2][4]≈M*Pset[2][4]*M' ? OK() : OH(name)

    name="congruence (IV)"; newTest(name)
    Pset1=randP(4, 2); # generate 2 positive definite 4x4 matrices
    Pset2=randP(4, 2);
    Pset=ℍVector₂([Pset1, Pset2]);
    U=𝕄Vector([randU(4), randU(4)])
    Qset=cong(U, Pset, MatrixVector₂)
    Qset[1][1]≈U[1]*Pset[1][1]*U[1]' &&
    Qset[1][2]≈U[1]*Pset[1][2]*U[2]' &&
    Qset[2][1]≈U[2]*Pset[2][1]*U[1]' &&
    Qset[2][2]≈U[2]*Pset[2][2]*U[2]' ? OK() : OH(name)

    ## 8. Spectral decompositions of positive matrices

    name="evd"; newTest(name)
    (Λ, U) = evd(P)
    U*Λ*U'≈P ? OK() : OH(name*" Real Input")
    (Λ, U) = evd(PC)
    U*Λ*U'≈PC ? OK() : OH(name*" Complex Input")


    name="frf"; newTest(name)
    F = frf(P)
    F*F'≈P ? OK() : OH(name*" Real Input")
    F = frf(PC)
    F*F'≈PC ? OK() : OH(name*" Complex Input")

    name="invfrf"; newTest(name)
    F = invfrf(P)
    F*P*F'≈I ? OK() : OH(name*" Real Input")
    F = invfrf(PC)
    F*PC*F'≈I ? OK() : OH(name*" Complex Input")

    name="spectralFunctions"; newTest(name);
    spectralFunctions(P, x->x+1); RUN()
    spectralFunctions(PC, abs2); RUN()
    spectralFunctions(D_, x->x+1); RUN()


    name="pow"; newTest(name)
    P½, P½ⁱ=pow(P_, 0.5, -0.5)
    P½*P½ⁱ≈I && P½*P½≈P_ ?  OK() : OH(name*" Real Input method 1")
    P½, P½ⁱ=pow(PC_, 0.5, -0.5)
    P½*P½ⁱ≈I && P½*P½≈PC_ ? OK() : OH(name*" Complex Input method 1")
    D½, D½ⁱ=pow(D_, 0.5, -0.5)
    D½*D½ⁱ≈I && D½*D½≈D_ ?  OK() : OH(name*" Real Input method 2")


    name="invsqrt"; newTest(name)
    P½ⁱ=invsqrt(P_)
    P½ⁱ*P_*P½ⁱ'≈I ? OK() : OH(name*" Real Input method 1")
    P½ⁱ=invsqrt(PC_)
    P½ⁱ*PC_*P½ⁱ'≈I ? OK() : OH(name*" Complex Input method 1")
    D½ⁱ=invsqrt(D_)
    D½ⁱ*D_*D½ⁱ'≈I ? OK() : OH(name*" Real Input method 2")


    name="sqr"; newTest(name)
    P²=sqr(P_)
    P² ≈ P_*P_' ? OK() : OH(name*" Real Input method 1")
    P²=sqr(PC_)
    P² ≈ PC_*PC_' ? OK() : OH(name*" Complex Input method 1")
    D²=sqr(D_)
    D² ≈ D_*D_' ? OK() : OH(name*" Real Input method 2")


    name="powerIterations"; newTest(name)
    Λ, U, iterations, convergence=powIter(P_, size(P_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 1")
    Λ, U, iterations, convergence=powIter(PC_, size(PC_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(PC_) && U*Λ*U'≈PC_ ? OK() : OH(name*" Complex Input Method 1")
    Λ, U, iterations, convergence=powIter(ℍ(P_), size(P_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 2")
    Λ, U, iterations, convergence=powIter(ℍ(PC_), size(PC_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(PC_) && U*Λ*U'≈PC_ ? OK() : OH(name*" Complex Input Method 2")
    Λ, U, iterations, convergence=powIter(LowerTriangular(P_), size(P_, 2), evalues=true)
    sort(diag(Λ))≈eigvals(P_) && U*Λ*U'≈P_ ? OK() : OH(name*" Real Input Method 3")


    name="choL"; newTest(name)
    L=choL(P_)
    L*L'≈P_ ? OK() : OH(name*" Real Input method 1")
    L=choL(PC_)
    L*L'≈PC_ ? OK() : OH(name*" Complex Input method 2")
    L=choL(D_)
    L*L'≈D_ ? OK() : OH(name*" Real Input method 2")


    name="choInv"; newTest(name)
    L, Linv=choInv(P_)
    invP_=inv(P_)
    L*L'≈P_ ? OK() : OH(name*" Real Input test L")
    Linv*Linv'≈invP_ ? OK() : OH(name*" Real Input test Linv")
    L, D, Linv=choInv(P_; kind=:LDLt)
    L*D*L'≈P_ ? OK() : OH(name*" Real Input, kind :LDLt test L")
    Linv*inv(D)*Linv'≈invP_ ? OK() : OH(name*" Real Input, kind :LDLt test Linv")

    L, Linv=choInv(PC_)
    invPC_=inv(PC_)
    L*L'≈PC_ ? OK() : OH(name*" Complex Input test L")
    Linv*Linv'≈invPC_ ? OK() : OH(name*" Complex Input test Linv")
    L, D, Linv=choInv(PC_; kind=:LDLt)
    L*D*L'≈PC_ ? OK() : OH(name*" Complex Input, kind :LDLt test L")
    Linv*inv(D)*Linv'≈invPC_ ? OK() : OH(name*" Complex Input, kind :LDLt test Linv")


    # 9. Decompositions involving triangular matrices

    # functions in SignalProcessinginP.jl
    println(" ")
    print("\n- Unit 'signalProcessing.jl'")

    name="randλ"; newTest(name);
    randλ(10); RUN()



    name="randΛ"; newTest(name);
    randΛ(10); RUN()
    randΛ(10, 2); RUN()



    name="randU"; newTest(name);
    U=randU(10)
    U'*U≈I ? OK() : OH(name*" Real Input")
    U=randU(ComplexF64, 10)
    U'*U≈I ? OK() : OH(name*" Complex Input")


    name="randP"; newTest(name);
    randP(3); RUN()
    randP(ComplexF64, 3); RUN()


    name="regularize!"; newTest(name)
    signalVar=tr(P)
    regularize!(P, SNR=10)
    signalPlusNoiseVar=tr(P)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 1")

    signalVar=tr(PC)
    regularize!(PC, SNR=10)
    signalPlusNoiseVar=tr(PC)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Complex Input Method 1")

    𝐏2=randP(5, 20)
    signalVar=sum(tr(P) for P in 𝐏2)
    regularize!(𝐏2, SNR=10)
    signalPlusNoiseVar=sum(tr(P) for P in 𝐏2)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")

    𝐏C2=randP(ComplexF64, 5, 20)
    signalVar=sum(tr(P) for P in 𝐏C2)
    regularize!(𝐏C2, SNR=10)
    signalPlusNoiseVar=sum(tr(P) for P in 𝐏C2)
    signalVar/(signalPlusNoiseVar-signalVar) ≈ 10 ? OK() : OH(name*" Real Input Method 2")


    name="gram"; newTest(name);
    gram(T); RUN()
    gram(W); RUN()


    name="trade"; newTest(name);
    trade(P); RUN()
    trade(PC); RUN()


    # functions in RiemannianGeometry.jl
    println(" ")
    print("\n- Unit 'riemannianGeometry.jl'")


    name="geodesic"; newTest(name);
    (geodesic(m, P, Q, 0.5) for m in metrics if m≠9); RUN()
    (geodesic(m, PC, QC, 0.5) for m in metrics if m≠9); RUN()
    (geodesic(m, 𝐃[1], 𝐃[2], 0.5) for m in metrics if m≠9); RUN()


    name="distanceSqr (I)"; newTest(name);
    for m in metrics distanceSqr(m, P) end; RUN()
    for m in metrics distanceSqr(m, P, Q) end; RUN()
    for m in metrics distanceSqr(m, PC) end; RUN()
    for m in metrics distanceSqr(m, PC, QC) end; RUN()

    name="distanceSqr (II)"; newTest(name);
    for m in metrics
         distanceSqr(m, D_)≈distanceSqr(m, ℍ(𝕄(D_))) ? OK() : OH(name*", metric "*string(m))
    end

    name="distanceSqr (III)"; newTest(name);
    for m in metrics
         distanceSqr(m, D_, E_)≈distanceSqr(m, ℍ(𝕄(D_)), ℍ(𝕄(E_))) ? OK() : OH(name*", metric "*string(m))
    end


    name="distance (I)"; newTest(name);
    # since it calls distanceSqr just check the call works
    for m in metrics distance(m, P) end; RUN()
    for m in metrics distance(m, P, Q) end; RUN()
    for m in metrics distance(m, PC) end; RUN()
    for m in metrics distance(m, PC, QC); end; RUN()

    name="distance (II)"; newTest(name);
    for m in metrics distance(m, D_) end; RUN()
    for m in metrics distance(m, D_, E_) end; RUN()

    𝐏=randP(8, 4)
    𝐏C=randP(ComplexF64, 8, 4)
    k=length(𝐏)
    kC=length(𝐏C)

    name="distanceSqrMat (I)"; newTest(name);
    for m in metrics
            L1=distanceSqrMat(Float64, m, 𝐏, ⏩=false)
            manualL1=𝕃{Float64}(diagm(0 => zeros(k)))
            for j=1:k-1, i=j+1:k manualL1[i, j]=distanceSqr(m, 𝐏[i], 𝐏[j]) end
            manualL1≈L1 ? OK() : OH(name*" Real Input, metric "*string(m))
    end

    name="distanceSqrMat (I ⏩ )"; newTest(name);
    for m in metrics
            L2=distanceSqrMat(Float64, m, 𝐏)
            manualL2=𝕃{Float64}(diagm(0 => zeros(k)))
            for j=1:k-1, i=j+1:k manualL2[i, j]=distanceSqr(m, 𝐏[i], 𝐏[j]) end
            manualL2≈L2 ? OK() : OH(name*" Real Input, metric "*string(m))
    end

    name="distanceSqrMat (II)"; newTest(name);
    for m in metrics
            L3=distanceSqrMat(Float64, m, 𝐏C, ⏩=false)
            manualL3=𝕃{Float64}(diagm(0 => zeros(kC)))
            for j=1:kC-1, i=j+1:kC manualL3[i, j]=distanceSqr(m, 𝐏C[i], 𝐏C[j]) end
            manualL3≈L3 ? OK() : OH(name*" Complex Input, metric "*string(m))
    end

    name="distanceSqrMat (II ⏩ )"; newTest(name);
    for m in metrics
            L4=distanceSqrMat(Float64, m, 𝐏C)
            manualL4=𝕃{Float64}(diagm(0 => zeros(kC)))
            for j=1:kC-1, i=j+1:kC manualL4[i, j]=distanceSqr(m, 𝐏C[i], 𝐏C[j]) end
            manualL4≈L4 ? OK() : OH(name*" Complex Input, metric "*string(m))
    end


    name="distanceMat (I)"; newTest(name);
    for m in metrics
            L5=distanceMat(Float64, m, 𝐏, ⏩=false)
            manualL5=𝕃{Float64}(diagm(0 => zeros(k)))
            for j=1:k-1, i=j+1:k manualL5[i, j]=distance(m, 𝐏[i], 𝐏[j]) end
            manualL5≈L5 ? OK() : OH(name*" Real Input, metric "*string(m))
    end

    name="distanceMat (I ⏩ )"; newTest(name);
    for m in metrics
            L6=distanceMat(Float64, m, 𝐏)
            manualL6=𝕃{Float64}(diagm(0 => zeros(k)))
            for j=1:k-1, i=j+1:k manualL6[i, j]=distance(m, 𝐏[i], 𝐏[j]) end
            manualL6≈L6 ? OK() : OH(name*" Real Input, metric "*string(m))
    end

    name="distanceMat (II)"; newTest(name);
    for m in metrics
            L7=distanceMat(Float64, m, 𝐏C, ⏩=false)
            manualL7=𝕃{Float64}(diagm(0 => zeros(kC)))
            for j=1:kC-1, i=j+1:kC manualL7[i, j]=distance(m, 𝐏C[i], 𝐏C[j]) end
            manualL7≈L7 ? OK() : OH(name*" Complex Input, metric "*string(m))
    end

    name="distanceMat (II ⏩ )"; newTest(name);
    for m in metrics
            L8=distanceMat(Float64, m, 𝐏C)
            manualL8=𝕃{Float64}(diagm(0 => zeros(kC)))
            for j=1:kC-1, i=j+1:kC manualL8[i, j]=distance(m, 𝐏C[i], 𝐏C[j]) end
            manualL8≈L8 ? OK() : OH(name*" Complex Input, metric "*string(m))
    end



    name="laplacian"; newTest(name);
    Dsqr=distanceSqrMat(logEuclidean, 𝐏)
    lap=laplacian(Dsqr); RUN()


    name="laplacianEigenMaps"; newTest(name);
    laplacianEM(lap, 2); RUN()


    name="spectralEmbedding"; newTest(name);
    spectralEmbedding(logEuclidean, 𝐏, 2); RUN()


    name="mean (I)"; newTest(name);
    for m=1:length(metrics)
            if m ∉ (7, 9) mean(metrics[m], P, Q) end end; RUN()
    for m=1:length(metrics)
            if m ∉ (6, 7, 9, 10) mean(metrics[m], 𝐏) end end; RUN()
    for m=1:length(metrics)
            if m ∉ (7, 9) mean(metrics[m], PC, QC) end end; RUN()
    for m=1:length(metrics)
            if m ∉ (6, 7, 9, 10) mean(metrics[m], 𝐏C) end end; RUN()
    for m=1:length(metrics)
            if m ∉ (6, 7, 9, 10) mean(metrics[m], 𝐃) end end; RUN()

    name="mean (II)"; newTest(name);
    k=length(𝐃)
    for m=1:length(metrics)
        if m ∉ (6, 7, 9, 10)
            D1=mean(metrics[m], 𝐃, ⏩=false)
            𝐃H=Vector{Hermitian}(undef, k)
            for i=1:k 𝐃H[i]=Hermitian(Matrix(𝐃[i])) end
            D2=mean(metrics[m], 𝐃H, ⏩=false)
            norm(𝕄(D1)-𝕄(D2))/k<0.0001 ? OK() : OH(name*" Real Diagonal Input, metric "*string(m))
        end
    end

    name="mean (⏩ )"; newTest(name);
    for m=1:length(metrics)
            if m ∉ (6, 7, 9, 10) mean(metrics[m], 𝐏) end end; RUN()
    for m=1:length(metrics)
            if m ∉ (7, 9) mean(metrics[m], 𝐃) end end; RUN()


    name="means"; newTest(name);
    means(logEuclidean, ℍVector₂([𝐏, 𝐐]); ⏩=false); RUN()
    means(logEuclidean, ℍVector₂([𝐏C, 𝐐C]); ⏩=false); RUN()
    means(logEuclidean, 𝔻Vector₂([𝐃, 𝐄]); ⏩=false); RUN()

    name="means (⏩ )"; newTest(name);
    means(logEuclidean, ℍVector₂([𝐏, 𝐐])); RUN()
    means(logEuclidean, ℍVector₂([𝐏C, 𝐐C])); RUN()
    means(logEuclidean, 𝔻Vector₂([𝐃, 𝐄])); RUN()


    name="generalizedMean"; newTest(name);
    𝐏2=ℍVector([P_, Q_])
    w=[0.2, 0.8]
    p=0.5
    ℍ( (P_^p+Q_^p)/2) ^(1/p) ≈ generalizedMean(𝐏2, p) ? OK() : OH(name*" Real Input 1")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w) ? OK() : OH(name*" Real Input 3")
    ℍ( (ℍ(0.2*P_^p)+ℍ(0.8*Q_^p))  )^(1/p) ≉ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 4")
    ℍ( (ℍ(0.4*P_^p)+ℍ(1.6*Q_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Real Input 5")
    𝐏2=ℍVector([PC_, QC_])
    w=[0.2, 0.8]
    ℍ( (PC_^p+QC_^p)/2) ^(1/p) ≈ generalizedMean(𝐏2, p) ? OK() : OH(name*" Complex Input 1")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 2")
    w=w.*2.0
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w) ? OK() : OH(name*" Complex Input 3")
    ℍ( (ℍ(0.2*PC_^p)+ℍ(0.8*QC_^p))  )^(1/p) ≉ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 4")
    ℍ( (ℍ(0.4*PC_^p)+ℍ(1.6*QC_^p))  )^(1/p) ≈ generalizedMean(𝐏2, p; w=w, ✓w=false) ? OK() : OH(name*" Complex Input 5")
    ((𝐃[1]^p+𝐃[2]^p)/2)^(1/p) ≈ generalizedMean(𝔻Vector([𝐃[1], 𝐃[2]]), p) ? OK() : OH(name*" Real Diagonal Input")

    name="generalizedMean"; newTest(name);
    generalizedMean(𝐏, 0.5; ⏩=false); RUN()
    generalizedMean(𝐏, 0.5; w=weights, ✓w=false, ⏩=false); RUN()
    generalizedMean(𝐏C, 0.5; ⏩=false); RUN()
    generalizedMean(𝐏C, 0.5; w=weights, ✓w=false, ⏩=false); RUN()
    generalizedMean(𝐃, 0.5; ⏩=false); RUN()
    generalizedMean(𝐃, 0.5; w=weights, ✓w=false, ⏩=false); RUN()


    name="geometricMean(⏩ )"; newTest(name);
    geometricMean(𝐏); RUN()
    geometricMean(𝐏, w=weights, ✓w=false);  RUN()
    geometricMean(𝐏C); RUN()
    geometricMean(𝐏C, w=weights, ✓w=false); RUN()
    geometricMean(𝐃); RUN()
    geometricMean(𝐃, w=weights, ✓w=false); RUN()

    name="geometricMean"; newTest(name);
    geometricMean(𝐏; ⏩=false); RUN()
    geometricMean(𝐏; w=weights, ✓w=false, ⏩=false); RUN()
    geometricMean(𝐏C; ⏩=false); RUN()
    geometricMean(𝐏C; w=weights, ✓w=false, ⏩=false); RUN()
    geometricMean(𝐃; ⏩=false); RUN()
    geometricMean(𝐃, w=weights, ✓w=false, ⏩=false); RUN()

    name="logdet0Mean"; newTest(name);
    w=[0.5, 0.5]
    P½, P½ⁱ=pow(P_, 0.5, -0.5)
    GM=P½*(P½ⁱ*Q_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(ℍVector([P_, Q_]); ⏩=false) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 1")
    ldG, iter, conv = logdet0Mean(ℍVector([P_, Q_]); w=w, ⏩=false) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Real Input 2")
    P½, P½ⁱ=pow(PC_, 0.5, -0.5)
    GM=P½*(P½ⁱ*QC_*P½ⁱ)^0.5*P½  # Fisher mean for k=2
    ldG, iter, conv = logdet0Mean(ℍVector([PC_, QC_]); ⏩=false) # logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 1")
    ldG, iter, conv = logdet0Mean(ℍVector([PC_, QC_]); w=w, ⏩=false) # weighted logdet0 mean for k=2
    GM ≈ ldG ? OK() : OH(name*" Complex Input 2")
    GM=exp(0.5*(log(𝐃[1])+log(𝐃[2])))  # Fisher mean for k=2 (=lofdet0 mean)
    ldG, iter, conv = logdet0Mean(𝔻Vector([𝐃[1], 𝐃[2]]); ⏩=false) # logdet0 mean for k=2
    isapprox(GM, ldG; atol=0.0001) ? OK() : OH(name*" Real Diagonal Input")

    name="logdet0Mean(⏩ )"; newTest(name);
    logdet0Mean(𝐏); RUN()
    logdet0Mean(𝐏; w=weights, ✓w=false); RUN()
    logdet0Mean(𝐏C); RUN()
    logdet0Mean(𝐏C; w=weights, ✓w=false); RUN()
    logdet0Mean(𝐃); RUN()
    logdet0Mean(𝐃, w=weights, ✓w=false); RUN()


    name="wasMean(⏩ )"; newTest(name);
    wasMean(𝐏); RUN()
    wasMean(𝐏; w=weights); RUN()
    wasMean(𝐏C); RUN()
    wasMean(𝐏C; w=weights); RUN()
    wasMean(𝐃); RUN()
    wasMean(𝐃; w=weights); RUN()

    name="wasMean"; newTest(name);
    wasMean(𝐏; ⏩=false); RUN()
    wasMean(𝐏; w=weights, ✓w=false, ⏩=false); RUN()
    wasMean(𝐏C; ⏩=false); RUN()
    wasMean(𝐏C; w=weights, ✓w=false, ⏩=false); RUN()
    wasMean(𝐃; ⏩=false); RUN()
    wasMean(𝐃, w=weights, ✓w=false, ⏩=false); RUN()



    name="powerMean(⏩ )"; newTest(name);
    powerMean(𝐏, 0.5); RUN()
    powerMean(𝐏, 0.5; w=weights); RUN()
    powerMean(𝐏C, 0.5); RUN()
    powerMean(𝐏C, 0.5; w=weights); RUN()
    powerMean(𝐃, 0.5); RUN()
    powerMean(𝐃, 0.5; w=weights); RUN()

    name="powerMean(⏩ )"; newTest(name);
    powerMean(𝐏, 0.5; ⏩=false); RUN()
    powerMean(𝐏, 0.5; w=weights, ✓w=false, ⏩=false); RUN()
    powerMean(𝐏C, 0.5; ⏩=false); RUN()
    powerMean(𝐏C, 0.5; w=weights, ✓w=false, ⏩=false); RUN()
    powerMean(𝐃, 0.5; ⏩=false); RUN()
    powerMean(𝐃, 0.5, w=weights, ✓w=false, ⏩=false); RUN()


    name="logMap"; newTest(name);
    logMap(Fisher, P, Q); RUN()
    logMap(Fisher, PC, QC); RUN()


    name="expMap"; newTest(name);
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


    name="vecP"; newTest(name);
    v=vecP(P); RUN()
    vC=vecP(PC); RUN()


    name="matP"; newTest(name);
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

    name="procrustes"; newTest(name);
    U=randU(size(P, 1))
    rotatedP=Hermitian(U'*P*U)
    procrustesU=procrustes(P, rotatedP)
    newP=Hermitian(procrustesU'*rotatedP*procrustesU)
    newP ≈ P ? OK() : OH(name*" Real Input")

    U=randU(ComplexF64, size(PC, 1))
    rotatedP=Hermitian(U'*PC*U)
    procrustesU=procrustes(PC, rotatedP)
    newP=Hermitian(procrustesU'*rotatedP*procrustesU)
    newP ≈ PC ? OK() : OH(name*" Complex Input")


    name="phi"; newTest(name);
    P=randP(4);
    Q=randP(4);
    geo=phi(P, Q, x->x^0.1; fispos=true);
    g=geodesic(Fisher, P, Q, 0.1);
    g≈geo ? OK() : OH(name*" Real Input");

    logmap=phi(P, Q, log);
    l=logMap(Fisher, Q, P);
    l≈logmap ? OK() : OH(name*" Real Input 2");

    cP=randP(ComplexF64, 4);
    cQ=randP(ComplexF64, 4);
    geo=phi(cP, cQ, x->x^0.1; fispos=true);
    g=geodesic(Fisher, cP, cQ, 0.1);
    g≈geo ? OK() : OH(name*" Complex Input");

    logmap=phi(cP, cQ, log);
    l=logMap(Fisher, cQ, cP);
    l≈logmap ? OK() : OH(name*" Complex Input 2");




    # functions in classification.jl
    println(" ")
    print("\n- Unit 'statistics.jl'")

    name="softmax"; newTest(name);
    g=[1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0]
    h= [0.0236405430216, 0.0642616585105, 0.17468129856, 0.47483299974438,
        0.0236405430216, 0.0642616585105, 0.17468129856]
    hh=softmax(g)
    hh ≈ h ? OK() : OH(name)

    name="mean (scalar version)"; newTest(name);
    g=[1.0, 2.0, 3.0]
    for m=1:length(metrics)
            if m ∉ (7, 9) mean(metrics[m], g) end end; RUN()

    name="std (scalar version)"; newTest(name);
    g=[1.0, 2.0, 3.0]
    for m=1:length(metrics)
            if m in (Fisher, Euclidean) std(metrics[m], g) end end; RUN()


end # function tests

function newTest(name::String)
    sleep(0.025)
    println(" ")
    print(rpad(name*":", 24))
end

OK()=print("⭐ ")

RUN()=print("☆ ")

function OH(name::String)
    print("⛔ ")
    push!(failing_tests, name)
end

SKIP()=print(" skypped")

failing_tests=[]

function testall()
    println("\n⭐","  PosDefManifold testing utility ", "⭐\n")
    println("Starting tests...\n")
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
