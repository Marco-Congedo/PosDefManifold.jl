makedocs(sitename="PosDefManifold", pages =
[ "index.md", "introToRiemannianGeometry.md", "MainModule.md", "riemannianGeometry.md",
 "linearAlgebra.md", "signalProcessing.md", "test.md"])

 deploydocs(
    repo = "github.com/Marco-Congedo/PosDefManifold.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
