using Documenter

makedocs(sitename="PosDefManifold", pages =
[ "index.md", "introToRiemannianGeometry.md", "MainModule.md", "riemannianGeometry.md",
 "linearAlgebra.md", "signalProcessing.md", "test.md"])

 deploydocs(
    repo = "github.com/Marco-Congedo/PosDefManifold.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
