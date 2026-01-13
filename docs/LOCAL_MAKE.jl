# Build the documenattion locally
# This si useful to check the documentation before deploying

push!(LOAD_PATH,"../src/")
push!(LOAD_PATH, @__DIR__)
using Documenter, DocumenterTools, DocumenterCitations, DocumenterInterLinks
using PosDefManifold

makedocs(
   sitename="PosDefManifold",
   remotes = nothing, 
   format = Documenter.HTML(
        prettyurls = false,
        edit_link = nothing
    ), # ELIMINATE prettyurls for deploying
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

