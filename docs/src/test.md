# test.jl

Most functions in **PosDefManifold** are tested, both for real and complex data input. This unit declares the function `testall()` that performs all tests.

Some functions are fully tested, the others are just executed.
Unce you ran it, for each method of each function,
a ⭐ sign is printed if the test is succesful, while
a ⛔ sign is printed if the test is not succesful.
A ☆ sign is printed if the function has been executed correctly.

Tests on functions for which a multi-threated version exist are indicated by symbol ( ⏩ ).

If there are fails, the concerned functions will be listed as *warnings*.

Note that the first time you execute the test it will take some time as the code will be compiled.

This here below is the output of the `testall()` function
(v0.5.2) run on the 11th of May 2025:

⭐ PosDefManifold testing utility⭐

Starting tests...

- Unit 'linearAlgebra.jl' 
typeofMatrix:           ☆  
dim:                    ☆  
det1:                   ⭐ ⭐  
tr1:                    ⭐ ⭐  
nearestPosDef:          ☆ ☆  
normalizeCol!:          ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
ispos:                  ⭐ ⭐  
colProd:                ⭐ ⭐ ⭐ ⭐  
colNorm:                ⭐ ⭐  
sumOfSqr:               ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
sumOfSqrDiag:           ⭐ ⭐ ⭐  
sumOfSqrTril:           ⭐ ⭐  
tr:                     ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
quadraticForm:          ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
fidelity:               ☆ ☆  
fDiag:                  ⭐  
DiagOfProd:             ⭐ ⭐  
mgs:                    ⭐ ⭐  
fVec:                   ⭐ ⭐ ⭐ ⭐  
congruence (I):         ⭐  
congruence (II):        ⭐  
congruence (III):       ⭐  
congruence (IV):        ⭐  
evd:                    ⭐ ⭐  
frf:                    ⭐ ⭐  
invfrf:                 ⭐ ⭐  
spectralFunctions:      ☆ ☆ ☆  
pow:                    ⭐ ⭐ ⭐  
invsqrt:                ⭐ ⭐ ⭐  
sqr:                    ⭐ ⭐ ⭐  
powerIterations:        ⭐ ⭐ ⭐ ⭐ ⭐  
choL:                   ⭐ ⭐ ⭐  
choInv:                 ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  

- Unit 'signalProcessing.jl' 
randλ:                  ☆  
randΛ:                  ☆ ☆  
randU:                  ⭐ ⭐  
randP:                  ☆ ☆  
regularize!:            ⭐ ⭐ ⭐ ⭐  
gram:                   ☆ ☆  
trade:                  ☆ ☆  

- Unit 'riemannianGeometry.jl' 
geodesic:               ☆ ☆ ☆  
distanceSqr (I):        ☆ ☆ ☆ ☆  
distanceSqr (II):       ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceSqr (III):      ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distance (I):           ☆ ☆ ☆ ☆  
distance (II):          ☆ ☆  
distanceSqrMat (I):     ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceSqrMat (I ⏩ ): ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceSqrMat (II):    ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceSqrMat (II ⏩ ):⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceMat (I):        ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceMat (I ⏩ ):    ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceMat (II):       ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
distanceMat (II ⏩ ):   ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
laplacian:              ☆  
laplacianEigenMaps:     ☆  
spectralEmbedding:      ☆  
mean (I):               ☆ ☆ ☆ ☆ ☆  
mean (II):              ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
mean (⏩ ):             ☆ ☆  
means:                  ☆ ☆ ☆  
means (⏩ ):            ☆ ☆ ☆  
generalizedMean:        ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐ ⭐  
generalizedMean:        ☆ ☆ ☆ ☆ ☆ ☆  
geometricMean(⏩ ):     ☆ ☆ ☆ ☆ ☆ ☆  
geometricMean:          ☆ ☆ ☆ ☆ ☆ ☆  
logdet0Mean:            ⭐ ⭐ ⭐ ⭐ ⭐  
logdet0Mean(⏩ ):       ☆ ☆ ☆ ☆ ☆ ☆  
wasMean(⏩ ):           ☆ ☆ ☆ ☆ ☆ ☆  
wasMean:                ☆ ☆ ☆ ☆ ☆ ☆  
powerMean(⏩ ):         ☆ ☆ ☆ ☆ ☆ ☆  
powerMean(⏩ ):         ☆ ☆ ☆ ☆ ☆ ☆  
logMap:                 ☆ ☆  
expMap:                 ☆ ☆  
vecP:                   ☆ ☆  
matP:                   ☆ ☆  
procrustes:             ⭐ ⭐  

- Unit 'statistics.jl' 
softmax:                ⭐  
mean (scalar version):  ☆  
std (scalar version):   ☆

Info: All tests were succesful!