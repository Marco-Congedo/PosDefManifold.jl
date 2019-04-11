# 1) hit SHIFT+ENTER
# 2) in the REPL: CTRL+v, ENTER
begin
    projectName="PosDefManifold"
    juliaCodeDir=homedir()*"\\Documents\\Code\\julia\\"
    # this file will be copied once to become the logo
    logoPngFile=homedir()*"\\Pictures\\cenerentola.png"

    projectDir       =   juliaCodeDir*projectName*"\\"
    srcDir           =   projectDir*"src\\"
    MainModule       =   srcDir*projectName*".jl"
    docsDir          =   projectDir*"docs\\"
    docsSrcDir       =   docsDir*"src\\"
    docsSrcAssetsDir =   docsSrcDir*"assets\\"
    docsBuildDir     =   docsDir*"build\\"

    # creates directories
    ispath(srcDir)           || mkpath(srcDir)
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

    using Documenter, Statistics, Revise, PosDefManifold
    #using Revise
    #include(MainModule)
    cd(docsDir)

    #2) copy and run the following line from the REPL
    #chmod(projectDocDir, recursive=true)
    #chmod(projectDocBuildDir, recursive=true)
    if Base.HOME_PROJECT[] !== nothing
        Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
    end


    clipboard("""makedocs(sitename="PosDefManifold", pages = [ "index.md", "introToRiemannianGeometry.md", "MainModule.md", "riemannianGeometry.md", "linearAlgebra.md", "signalProcessing.md", "test.md"])""")
    @info("\nhit CTRL+V+ENTER on the REPL for building the documentation.");
end





#makedocs(sitename="PosDefManifold")

#makedocs(sitename="PosDefManifold", assets = ["assets/style.css"])

#=

makedocs(
 modules = [GeophysicalFlows],
 doctest = false,
   clean = truemakedocs(sitename="PosDefManifold



checkdocs = :all,
  format = :html,
 authors = "Gregory L. Wagner and Navid C. Constantinou",
sitename = "GeophysicalFlows.jl",
   pages = Any[
            "Home" => "index.md",
            "Modules" => Any[
              "modules/twodturb.md",
              "modules/barotropicqg.md",
              "modules/multilayerqg.md"
            ],
            "DocStrings" => Any[
            "man/types.md",
            "man/functions.md"]
           ]
)

deploydocs(repo = "github.com/FourierFlows/GeophysicalFlows.jl.git")


istravis = "TRAVIS" ∈ keys(ENV)

makedocs(
  format = Documenter.HTML(prettyurls=istravis),
  sitename = "GeoStats.jl",
  authors = "Júlio Hoffimann Mendes",
  assets = ["assets/style.css"],
  pages = [
    "Home" => "index.md",
    "User guide" => [
      "Problems & solvers" => "problems_and_solvers.md",
      "Spatial data" => "spatialdata.md",
      "Domains" => "domains.md",
      "Variography" => [
        "empirical_variograms.md",
        "theoretical_variograms.md",
        "fitting_variograms.md"
      ],
      "Kriging estimators" => "estimators.md",
      "Solver comparisons" => "comparisons.md",
      "Spatial statistics" => "statistics.md",
      "Plotting" => "plotting.md"
    ],
    "Tutorials" => "tutorials.md",
    "Contributing" => "contributing.md",
    "About" => [
      "Community" => "about/community.md",
      "License" => "about/license.md",
      "Citing" => "about/citing.md"
    ],
    "Developer guide" => "developers.md",
    "Index" => "links.md"
  ]
)

deploydocs(repo="github.com/juliohm/GeoStats.jl.git")
=#
