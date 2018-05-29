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
using FFD
using FactCheck

comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)


resize!(ARGS, 1)
ARGS[1] = "./input_vals_3d_parallel_np3and4.jl"

opts = PDESolver.read_input(ARGS[1])
sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)
orig_vert_coords = deepcopy(mesh.vert_coords)
geom_faces = opts["BC2"]

facts("--- Checking evaldXdControlPointProduct for 3D DG Mesh ---") do

  # for i = 1:comm_size
  #   if my_rank == i-1
  #     println("mesh.vert_sharing.vert_nums = \n$(mesh.vert_sharing.vert_nums)")
  #     println("mesh.vert_sharing.rev_mapping = \n$(mesh.vert_sharing.rev_mapping)\n")
  #   end
  #   MPI.Barrier(comm)
  # end

  # for i = 1:comm_size
  #   if my_rank == i-1
  #     println("Rank $my_rank,\nmesh.peer_parts = \n$(mesh.peer_parts)")
  #     println("shared interfaces = \n$(mesh.shared_interfaces)\n")
  #   end
  #   MPI.Barrier(comm)
  # end

  ndim = mesh.dim
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [4,4,2]
  offset = [0.5, 0.5, 0.5] # No offset in the X & Y direction
  geom_faces = opts["BC2"]

  # Create a mapping object using nonlinear mapping
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

  for i = 1:comm_size
    fname = string("controlPoints_rank", my_rank)
    writeControlPointsVTS(map, fname)
    MPI.Barrier(comm)
  end

  MPI.Barrier(comm)
  fill!(map.work, 0.0)

  # FFD.bndryMPIRanks(mesh, geom_faces)


  # Create seed vector
  # - Get original wall coordinates
  orig_wallCoords = FFD.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
  for i = 1:comm_size
    if my_rank == i-1
      println("on rank $my_rank, orig_wallCoords = \n$(orig_wallCoords)")
    end
    MPI.Barrier(comm)
  end


  Xs_bar = ones(3, size(orig_wallCoords,2))
  evaldXdControlPointProduct(map, mesh, vec(Xs_bar))

  fname = "./testvalues/evaldXdControlPointProduct_tet8cube.dat"

  test_values = readdlm(fname)
  @fact length(map.work) --> length(test_values)
  for i = 1:length(map.work)
    err = abs.(test_values[i] - map.work[i])
    @fact err --> less_than(1e-14) "problem at index $i"
  end
  
  #=
  for i = 1:comm_size
    fname2 = string("./testvalues/evaldXdControlPointProduct_tet_cube_8el_rank",
                    my_rank, ".dat")
    f = open(fname2, "w")
    for j = 1:length(map.work)
      println(f, map.work[j])
    end
    close(f)
  end
  =#

end # End facts("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")

