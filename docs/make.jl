push!(LOAD_PATH,"../src/")
using Documenter, PosDefManifold

makedocs (
   sitename="PosDefManifold",
   modules=[PosDefManifold],
   pages =
   [
      "index.md",
      "introToRiemannianGeometry.md",
      "MainModule.md",
      "riemannianGeometry.md",
      "linearAlgebra.md",
      "signalProcessing.md",
      "test.md"
   ]
)

deploydocs(
    repo = "github.com/Marco-Congedo/PosDefManifold.jl.git"
)
