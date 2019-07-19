#   This unit is not necessary for the PosDefManifold package.
#   It allows to install the PosDefManifold package
#   and to build the documentation locally.
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the PosDefPackage from github or Julai registry,
#      uninstall it.
#   2) Change the `juliaCodeDir` path here below to the path
#           where the PosDefMaifold folder is located on your computer.
#   3) Under Linux, replace all '\\' with `/`
#   4) Put the cursor in this unit and hit SHIFT+CTRL+ENTER
#
#   Nota Bene: all you need is actually the 'push' line and
#   the 'using' line. You can safely delete the rest once
#   you have identified the 'srcDir' to be used in the push command.

begin
    # change the 'juliaCodeDir' path to the folder where your projects are
    juliaCodeDir= homedir()*"\\Documents\\Code\\julia\\"

    projectName = "PosDefManifold"
    srcDir      = juliaCodeDir*projectName*"\\src\\"
    docsDir     = juliaCodeDir*projectName*"\\docs\\"

    push!(LOAD_PATH, srcDir)

    using Documenter, DocumenterTools, Statistics,
    Base.Threads, LinearAlgebra, BenchmarkTools, Revise,
    PosDefManifold

    # add moduls to be used
    Modules = juliaCodeDir*"Modules"
    push!(LOAD_PATH, Modules)

    using IOtxt

    # for compiling the documentation
    cd(docsDir)
    clipboard("""makedocs(sitename="$projectName", modules=[$projectName])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end
