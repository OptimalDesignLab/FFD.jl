# runtest.jl
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

using FFD  # have FFD initialize MPI
#using PDESolver
using PdePumiInterface
using SummationByParts
using ODLCommonTools
using ArrayViews
#using Utils
using MPI
# using MeshMovement
using Base.Test

# Include the actual tests
include("./test_b-splines.jl")
include("./test_serial_FFD_2D_mesh.jl")
include("test_2d.jl")
include("./test_contractWithdGdB.jl")
include("./test_constraints.jl")
#include("./test_serial_FFD_3D_mesh.jl")

# Debug

#MPI.Finalize()

Base.Test.exitstatus()
