#   This script is not part of the PosDefManifold package.
#   It allows to build the PosDefManifold package
#   and its documentation locally from the source code,
#   without actually installing the package.
#   It is used for developing purposes using the Julia
#   `Revise` package (that you need to have installed on your PC,
#   together with the `Documenter` package for building the documentation).
#   You won't need this script for using the package.
#
#   DIRECTIONS:
#   1) If you have installed the PosDefManifold package
#      from github or Julia registry, uninstall it.
#   2) Change the `juliaCodeDir` path here below to the path
#           where the PosDefManifold folder is located on your PC.
#   3) Under Linux, replace all '\\' with `/`
#   4) Put the cursor in this unit and hit SHIFT+CTRL+ENTER
#
#   Nota Bene: all you need for building the package is actually
#   the 'push' line and the 'using' line.
#   You can safely delete the rest once
#   you have identified the 'srcDir' to be used in the push command.

begin
    # change the 'juliaCodeDir' path to the folder where your projects are
    juliaCodeDir= homedir()*"\\Documents\\@ Documenti\\Code\\julia\\"
    projectName = "PosDefManifold"
    srcDir      = juliaCodeDir*projectName*"\\src\\"
    docsDir     = juliaCodeDir*projectName*"\\docs\\"

    push!(LOAD_PATH, srcDir)

    using Documenter, LinearAlgebra, Statistics, Base.Threads,
          Revise, PosDefManifold

    # add other local moduls to be used, e.g.,
    # Modules = juliaCodeDir*"Modules"
    # push!(LOAD_PATH, Modules)
    # using IOtxt
<<<<<<< HEAD

    # for compiling the documentation
    cd(docsDir);
    clipboard("""makedocs(sitename="PosDefManifold", modules=[PosDefManifold])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end

# for compiling the documentation
cd(docsDir);
clipboard("""makedocs(sitename="PosDefManifold", modules=[PosDefManifold])""")
@info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
#makedocs(sitename="PosDefManifold", modules=[PosDefManifold])

