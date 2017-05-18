# runtest.jl
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

using PdePumiInterface
using PDESolver
using SummationByParts
using ODLCommonTools
using ArrayViews
using Utils
using MPI
# using MeshMovement
using FreeFormDeformation
using FactCheck

# Include functions that will be used for performing tests
include("./test_functions.jl")

# Include the actual tests
include("./test_b-splines.jl")
include("./test_linearFFD.jl")
include("./test_FFD_derivatives.jl")
include("./test_constraints.jl")

# Debug

MPI.Finalize()
