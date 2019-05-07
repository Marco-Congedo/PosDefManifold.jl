# 1) hit SHIFT+ENTER
# 2) in the REPL: CTRL+v, ENTER

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

    # Copy logo png file
    path, ext=splitext(logoPngFile)
    destination=docsSrcAssetsDir*"logo.png"
    if isfile(logoPngFile) && !isfile(destination)
        if ext==".png" cp(logoPngFile, docsSrcAssetsDir*"logo.png")
        else @warn("The logo file is not a .png file.")
        end
    end

    #just
    wd=pwd()

    println("\n Your Machine: ",gethostname(),", in ",Sys.MACHINE)
    println(" Kernel: ",Sys.KERNEL,", with word size ",Sys.WORD_SIZE)
    println(" CPU Threads: ",Sys.CPU_THREADS, "\n")
    # Sys.BINDIR # julia bin directory

    #--------------------------------
    # for building the documentation:
    #--------------------------------

    #1) run the following three lines
    push!(LOAD_PATH, srcDir)

    using Documenter, DocumenterTools, Statistics,
    Revise,
    PosDefManifold, LinearAlgebra, BenchmarkTools
    cd(docsDir)

    #2) copy and run the following line from the REPL
    #chmod(projectDocDir, recursive=true)
    #chmod(projectDocBuildDir, recursive=true)
    if Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end

    clipboard("""makedocs(sitename="PosDefManifold", modules=[PosDefManifold])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end
