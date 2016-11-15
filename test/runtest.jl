# runtest.jl
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

using PdePumiInterface
using SummationByParts
using ODLCommonTools
using ArrayViews
using Utils
using MPI
using MeshMovement
using FreeFormDeformation
using FactCheck

include("./test_b-splines.jl")


#=
# Source file includes
push!(LOAD_PATH, "../src/")

using FactCheck
using FreeFormDeformation

# Create Test mesh for tests
nnodes = [3,3,3]  # Number of nodes of the FE grid that need to be mapped
nodes_xyz = zeros(nnodes[1], nnodes[2], nnodes[3], 3)
origin = [1,1,1]
incz = 0.0
for k = 1:nnodes[3]
  incy = 0.0
  for j = 1:nnodes[2]
    incx = 0.0
    for i = 1:nnodes[1]
      nodes_xyz[i,j,k,1] = origin[1] + incx
      nodes_xyz[i,j,k,2] = origin[2] + incy
      nodes_xyz[i,j,k,3] = origin[3] + incz
      incx += 1
    end
    incy += 1
  end
  incz += 1
end

# Create Mapping Object
ndim = 3
order = [2,2,2]  # Order of B-splines in the 3 directions
nControlPts = [3,3,3]
map = Mapping(ndim, order, nControlPts, nnodes)

# Create BoundingBox object
offset = [0.5,0.5,0.5]  # offset for the bounding box
geom_bounds = [1. 1. 1.;3. 3. 3.]
box = BoundingBox(ndim, geom_bounds, offset)

# Tests
include("./test_b-splines.jl")
include("./test_linearFFD.jl")
=#
