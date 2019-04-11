#   Unit Plots.jl, part of PosDefManifold Package for julia language
#   v 0.0.1 - last update 21 Mars 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DESCRIPTION
#   This Unit produces, show and save several plots to compare
#   the characteristics of the metrics implemented in PosDefManifold
#   to perform operations on the Riemannian manifold
#   of Symmetric Positive Definite (SPD) or Hermitian matrices
#
#   CONTENT
#   1. Comparing distances

#   Nota Bene: try calling plotly() to use this backend
#   if your default backend does not work


# BEGIN HERE
begin
 include("START HERE.jl")
 using Plots
 using LinearAlgebra
 using Statistics
end

function randQ(n::Int; df::Int=2, eigvalsSNR::Real=10e3)
     U=rand(n, n).-0.5
     return Symmetric(U * randEigvalsMat(n, df=df, eigvalsSNR=eigvalsSNR) * U')
     #return Symmetric(U * randP(n, df=df, eigvalsSNR=eigvalsSNR) * U')

end

#    1. Comparing distances


 # create and array of metrics of type Metric (declared in main module PosDefManifold.jl)
 # this is used for writing skinner code for both computations and for plots
 # The metrics are ordered in a particular way to highlight some results
 # The order is visible in the titles of the subplots (left to right and top to bottom)
order=(4, 5, 3, 6, 9, 2, 10, 8, 7, 1)
metrics=[Metric(i) for i in order]



 # set some variables
nm=length(metrics)  # # of metrics
 nk=100              # # of random points (matrices)
 dim=[8, 32, 128]    # comparisons will be done for three dimensions of matrices
 nn=length(dim)
 whiteNoise=10000         # regularization of generated random points expressed in SNR
 snr=100              # signal-to-noise ratio
df=2

 # generate random points and put them in 'P' and 'Q'
P=Array{Array}(undef, nn); # array [1:nn] of a nk,nm-dimensional array of matrices
 Q=Array{Array}(undef, nn); # array [1:nn] of a nk,nm-dimensional array of matrices
 for n=1:nn
  #R=randP(dim[n], nk, df=2, eigvalsSNR=whiteNoise, SNR=snr)
  R=[randQ(dim[n], df=df, eigvalsSNR=whiteNoise) for i=1:nk]  # generate a nk-dimensional array of random SPD matrices
  P[n]=R;
#  R=randP(dim[n], nk, df=2, eigvalsSNR=whiteNoise, SNR=snr)
  R=[randQ(dim[n], df=df, eigvalsSNR=whiteNoise) for i=1:nk]  # generate a nk-dimensional array of random SPD matrices
  Q[n]=R
 end

 ## Compute distances for all dimensions, all points and all metrics
X=[Rdistance(P[n][k], Q[n][k], metrics[m]) for n=1:nn, k=1:nk, m=1:nm]

# X=[Rnorm(P[n][k], metrics[m]) for n=1:nn, k=1:nk, m=1:nm]

 # Plot the Fisher distance on x-axis against all other distances (one for
 # each subplot) on the y-axis. A separate plot is done for each dimension

dispersions=[mean((X[n, :, 1]).^2)/dim[n] for n=1:nn]
for n=1:nn println(dispersions[n]) end

#spear=[cor(sortperm(X[n, :, 1]), sortperm(X[n, :, i])) for n=1:nn, i=1:nm]
spear=[cor(X[n, :, 1], X[n, :, i]) for n=1:nn, i=1:nm]


dig3(x)= round(x*100)/100

 # generate titles for subplots
 titles=[["$(metrics[i])"*" "*"$(dig3(spear[n, i]))" for j = 1:1, i=2:length(metrics)] for n=1:nn]

 for n=1:nn
    plots=[plot(X[n, :, 1], X[n, :, p+1],seriestype=:scatter) for p=1:nm-1]
    plot(plots..., title=titles[n], titlefontsize=10, legend=false, axis=nothing,
         markersize=3, markercolor=:red, markerstrokealpha=0.5, show=true,
         smooth=false)
    png(docsSrcAssetsDir*"Dist_n$(dim[n])"*"_k$nk"*"_wn$whiteNoise"*"_snr$snr")
 end

heatmap(Matrix(P[3][10]), yflip=true, c=:bluesreds)

#    2. Geodesics
#=
order=(4, 3, 2, 10, 8, 7, 1)
metrics=[Metric(i) for i in order]
metrics

# generate titles for subplots
titles=["$(metrics[i])" for j = 1:1, i=2:length(metrics)]

step=(0:0.1:1)


n=5
plot()
A=randn(n, n); P=Hermitian(A*A')
#P=randP(n, whiteNoise=1)
A=randn(n, n); Q=Hermitian(A*A')
#Q=inv(P)
#Q=randP(n, whiteNoise=1)
P
Q
 #𝔊=[Rmove(P, Q, a) for a in step]
𝔊=[Rmove(P, Q, a, metrics[5]) for a in step]

 #𝔊=[Rpow(P, a) for a in step]

gtr=[log(tr(G)/n) for G in 𝔊]

gdet=[logdet(G) for G in 𝔊]

plot!(gtr, gdet)

plot!(gtr, gdet, seriestype=:scatter)

gtr
gdet

P
Q

=#

#heatmap(P, yflip=true, c=:bluesreds)
#heatmap(P, yflip=true, c=:grays)
