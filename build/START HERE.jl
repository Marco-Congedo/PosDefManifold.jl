#   This unit allows to install the PosDefManifold package
#   and to build the documentation locally.
#   v 0.1.3 - last update 10th of Mai 2019
#
#   MIT License
#   Copyright (c) 2019, Marco Congedo, CNRS, Grenobe, France:
#   https://sites.google.com/site/marcocongedo/home
#
#   DIRECTIONS:
#   1) If you have installed the PosDefPackage from github, uninstall it.
#   2) Change the `juliaCodeDir` path here below to the path
#           where the PosDefMaifold folder is located on your computer.
#   3) Put the cursor in this unit and hit SHIFT+CTRL+ENTER

begin
    projectName="PosDefManifold"
    juliaCodeDir=homedir()*"\\Documents\\Code\\julia\\"
    # this file will be copied once to become the logo
    logoPngFile=homedir()*"\\Pictures\\cenerentola.png"

    projectDir       =   juliaCodeDir*projectName*"\\"
    srcDir           =   projectDir*"src\\"
    testDir          =   projectDir*"test\\"
    MainModule       =   srcDir*projectName*".jl"
    docsDir          =   projectDir*"docs\\"
    docsSrcDir       =   docsDir*"src\\"
    docsSrcAssetsDir =   docsSrcDir*"assets\\"
    docsBuildDir     =   docsDir*"build\\"

    # creates directories
    ispath(srcDir)           || mkpath(srcDir)
    ispath(testDir)          || mkpath(testDir)
    ispath(docsBuildDir)     || mkpath(docsBuildDir)
    ispath(docsSrcAssetsDir) || mkpath(docsSrcAssetsDir)

    # Copy logo png file if it does not exist
    path, ext=splitext(logoPngFile)
    destination=docsSrcAssetsDir*"logo.png"
    if isfile(logoPngFile) && !isfile(destination)
        if ext==".png" cp(logoPngFile, docsSrcAssetsDir*"logo.png")
        else @warn("The logo file is not a .png file.")
        end
    end

    # just an alias for convenience
    wd=pwd()

    #--------------------------------
    # for building the documentation:
    #--------------------------------

    push!(LOAD_PATH, srcDir)

    using Documenter, DocumenterTools, Statistics,
    Base.Threads, LinearAlgebra, BenchmarkTools, Revise,
    PosDefManifold

    cd(docsDir)

    #chmod(projectDocDir, recursive=true)
    #chmod(projectDocBuildDir, recursive=true)
    if  Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end

    clipboard("""makedocs(sitename="PosDefManifold", modules=[PosDefManifold])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end
