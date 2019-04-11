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
 using Plots, LinearAlgebra, Statistics
end


#    1. Comparing distances
begin

 # create and array of metrics of type Metric (declared in main module PosDefManifold.jl)
 # this is used for writing skinner code for both computations and for plots
 # The metrics are ordered in a particular way to highlight some results
 # The order is visible in the titles of the subplots (left to right and top to bottom)
 order=(4, 5, 3, 6, 9, 2, 10, 8, 7, 1)
 metrics=[Metric(i) for i in order]

 # generate titles for subplots
 titles=["$(metrics[i])" for j = 1:1, i=2:length(metrics)]

 # set some variables
 nm=length(metrics)  # # of metrics
 nk=100              # # of random points (matrices)
 dim=[2, 10, 100]    # comparisons will be done for three dimensions of matrices
 nn=length(dim)
 whiteNoise=10000          # regularization of generated random points expressed in SNR
 snr=0.1                  # signal-to-noise ratio
 df=2

 # generate random points and put them in 'P' and 'Q'
 P=Array{Array}(undef, nn); # array [1:nn] of a nk,nm-dimensional array of matrices
 Q=Array{Array}(undef, nn); # array [1:nn] of a nk,nm-dimensional array of matrices
 for n=1:nn
  R=randP(dim[n], nk, df=2, eigvalsSNR=whiteNoise, SNR=snr)  # generate a nk-dimensional array of random SPD matrices
  P[n]=R;
  R=randP(dim[n], nk, df=2, eigvalsSNR=whiteNoise, SNR=snr)
  Q[n]=R
 end

 ## Compute distances for all dimensions, all points and all metrics
 X=[Rdistance(P[n][k], Q[n][k], metrics[m]) for n=1:nn, k=1:nk, m=1:nm]


 # Plot the Fisher distance on x-axis against all other distances (one for
 # each subplot) on the y-axis. A separate plot is done for each dimension
 for n=1:nn
    plots=[plot(X[n, :, 1], X[n, :, p+1],seriestype=:scatter) for p=1:nm-1]
    plot(plots..., title=titles, titlefontsize=10, legend=false, axis=nothing,
         markersize=3, markercolor=:red, markerstrokealpha=0.5, show=true,
         smooth=true)
    png(docsSrcAssetsDir*"n$(dim[n])"*"_k$nk"*"_df$df"*"_wn$whiteNoise"*"_snr$snr")
 end
end

#    2. Geodesics
order=(4, 3, 2, 10, 8, 7, 1)
metrics=[Metric(i) for i in order]
metrics

# generate titles for subplots
titles=["$(metrics[i])" for j = 1:1, i=2:length(metrics)]

step=(0:0.1:1)

heatmap(Matrix(P[1]), yflip=true, c=:bluesreds)
#heatmap(P, yflip=true, c=:grays)
