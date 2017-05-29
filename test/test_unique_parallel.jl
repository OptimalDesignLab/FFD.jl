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
geom_faces = opts["BC2"]

wallCoords = FreeFormDeformation.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
local_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, geom_faces)

# for i = 1:comm_size
#   if my_rank == i-1
#     println("local unique wallCoords rank $(i-1) =\n$local_wallCoords")
#     println("wallCoords rank $(i-1) =\n$wallCoords")
#   end
#   MPI.Barrier(comm)
# end

facts("--- Check if Unique wall coordinates are being computed correctly across all ranks ---") do
  if my_rank == 0
    @fact wallCoords --> roughly([0.0 0.0 0.0 0.0 0.0 0.0
                                  0.0 0.5 0.0 0.5 1.0 1.0
                                  0.5 1.0 1.0 0.5 1.0 0.5], atol=1e-14)
  end

  if my_rank == 1
    @fact wallCoords --> roughly([0.0 0.0 0.0
                                  0.0 0.5 1.0
                                  0.0 0.0 0.0], atol=1e-14)
  end

end # End facts

facts("--- Checking evaldXdControlPointProduct for 3D DG Mesh ---") do

  ndim = mesh.dim
  order = [2,2,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,2,2]
  offset = [0.1, 0.0, 0.0] # No offset in the X & Y direction
  geom_faces = opts["BC2"]

  # Create a mapping object using nonlinear mapping
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

  for i = 1:comm_size
    fname = string("controlPoints_rank", my_rank)
    writeControlPointsVTS(map, fname)
  end

  fill!(map.work, 0.0)

  # Create seed vector
  # - Get original wall coordinates
  orig_wallCoords = FreeFormDeformation.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
  Xs_bar = ones(3, size(orig_wallCoords,2))
  evaldXdControlPointProduct(map, mesh, vec(Xs_bar))

  fname = "./testvalues/evaldXdControlPointProduct_tet8cube.dat"

  if my_rank == 0
    println("rank = $my_rank\n")
    test_values = readdlm(fname)
    @fact length(map.work) --> length(test_values)

    fname2 = "./testvalues/evaldXdControlPointProduct_tet8cube_rank0.dat"
    f = open(fname2, "w")
    for i = 1:length(map.work)
      println(f, map.work[i])
    end
    close(f)
    # for i = 1:length(map.work)
    #   err = abs(test_values[i] - map.work[i])
    #   @fact err --> less_than(1e-14) "problem at index $i"
    # end

  end

  MPI.Barrier(comm)

  if my_rank == 1

    println("\nrank = $my_rank\n")

    test_values = readdlm(fname)
    @fact length(map.work) --> length(test_values)

    fname2 = "./testvalues/evaldXdControlPointProduct_tet8cube_rank1.dat"
    f = open(fname2, "w")
    for i = 1:length(map.work)
      println(f, map.work[i])
    end
    close(f)
  end
end # End facts("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")

MPI.Finalize()
