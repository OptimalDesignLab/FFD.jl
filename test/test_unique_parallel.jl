# MPI Declarations
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

if !MPI.Initialized()
  MPI.Init()
end
comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)


resize!(ARGS, 1)
ARGS[1] = "./input_vals_3d_parallel.jl"

opts = PDESolver.read_input(ARGS[1])
sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)
orig_vert_coords = deepcopy(mesh.vert_coords)

MPI.Barrier(comm)
for i = 1:comm_size
  if my_rank == i-1
    println("VertSharing on rank $i")
    println("rank $i, mesh.VertSharing.npeers = $(mesh.vert_sharing.npeers)")
    println("rank $i, mesh.VertSharing.peer_nums = $(mesh.vert_sharing.peer_nums)")
    println("rank $i, mesh.VertSharing.counts = $(mesh.vert_sharing.counts)")
    println("rank $i, mesh.VertSharing.vert_nums = $(mesh.vert_sharing.vert_nums)")
    println("rank $i, mesh.VertSharing.rev_mapping = $(mesh.vert_sharing.rev_mapping)")
  end
end

ndim = mesh.dim
order = [4,4,2]  # Order of B-splines in the 3 directions
nControlPts = [4,4,2]
offset = [0.1, 0.1, 0.1] # No offset in the X & Y direction
geom_faces = opts["BC2"]

# Create a mapping object using nonlinear mapping
map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)
if my_rank == 0
  writeControlPointsVTS(map)
end

wallCoords = FreeFormDeformation.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
local_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, geom_faces)

for i = 1:comm_size
  if my_rank == i-1
    println("local unique wallCoords rank $(i-1) =\n$local_wallCoords")
    println("wallCoords rank $(i-1) =\n$wallCoords")
  end
  MPI.Barrier(comm)
end
