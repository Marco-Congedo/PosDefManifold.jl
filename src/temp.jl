using LinearAlgebra, PosDefManifold

function MDM(data, metrics=(Fisher, Fisher))
    if length(metrics)==1
        distmetric=metrics[1]
        meanmetric=metrics[1]
    else
        distmetric=metrics[1]
        meanmetric=metrics[2]
    end
    return 
end

n=2
k=2
t=(20, 30)
C1=randP(n)
C2=randP(n)

℘=ℍVector(undef, t[1])
for k in 1:t[1] ℘[k]=geodesic(randP(n), C1, 0.5) end

ℛ=ℍVector(undef, t[2])
for k in 1:t[2] ℛ[k]=geodesic(randP(n), C2, 0.5) end

means=MDM((℘, ℛ))

means[1]
means[2]
