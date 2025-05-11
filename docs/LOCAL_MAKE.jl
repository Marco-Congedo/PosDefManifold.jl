# Build the documenattion locally
# This si useful to check the documentation before deploying

push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, @__DIR__)
using Documenter, DocumenterTools, DocumenterCitations, DocumenterInterLinks
using PosDefManifold

makedocs(
   sitename="PosDefManifold",
   format = Documenter.HTML((prettyurls = false)), # ELIMINATE pretty URL for deploying
   authors="Marco Congedo, CNRS, Grenoble, France",
   modules=[PosDefManifold],
   pages =  [
      "index.md",
      "introToRiemannianGeometry.md",
      "MainModule.md",
      "riemannianGeometry.md",
      "linearAlgebra.md",
      "statistics.md",
      "signalProcessing.md",
      "test.md",
   ]
)

deploydocs(
    repo = "github.com/Marco-Congedo/PosDefManifold.jl.git",
    target = "build",
    devurl = "dev",
)
