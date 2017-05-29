# All tests for 3D parallel PUMI meshes go here

# MPI Declarations
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

facts("--- Check if Unique wall coordinates are being computed correctly across all ranks ---") do

  geom_faces = opts["BC2"]
  wallCoords = FreeFormDeformation.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
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
MPI.Barrier(comm)

facts("--- Checking FFD on 3D parallel DG Pumi meshes ---") do

  context("Check control point manipulation with nonlinear mapping on full mesh") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    full_geom = true

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, full_geom)
    @fact box.box_bounds --> roughly([0.0 0.0 -0.5
                                      1.0 1.0 1.5], atol=1e-14)

    # Rigid body rotation
    theta = -20*pi/180  # Rotate wall coordinates by 20 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
    # Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # Rigid body translation
    map.cp_xyz[1,:,:,:] += 0.2
    map.cp_xyz[2,:,:,:] += 0.3
    map.cp_xyz[3,:,:,:] += 0.5

    vertices = evalVolume(map, mesh)
    commitToPumi(map, mesh, sbp, vertices)

    fname = string("./testvalues/translation_plus_rotation_DG_3D_tet8cube_rank", my_rank,".dat")
    # f = open(fname, "w")
    # for i = 1:length(mesh.vert_coords)
    #   println(f, mesh.vert_coords[i])
    # end
    # close(f)

    test_values = readdlm(fname)
    for i = 1:length(test_values)
      err = abs(test_values[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

  end # End context("Check control point manipulation with nonlinear mapping on full mesh")

  MPI.Barrier(comm)
  # Reset the coordinates and mesh to the original value
  for i = 1:mesh.numEl
    update_coords(mesh, i, orig_vert_coords[:,:,i])
  end
  commit_coords(mesh, sbp)

  context("Check control point manipulation on a geometry face") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0.5, 0.5, 0.5] # No offset in the X & Y direction
    geom_faces = opts["BC2"]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

    # Rigid body rotation
    theta = -2*pi/180  # Rotate wall coordinates by 2 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
    # Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # Rigid body translation
    map.cp_xyz[1,:,:,:] += 0.02
    map.cp_xyz[2,:,:,:] += 0.03
    map.cp_xyz[3,:,:,:] += 0.05

    vertices = evalSurface(map, mesh)
    commitToPumi(map, mesh, sbp, vertices)

    fname = string("./testvalues/translation_plus_rotation_DG_3D_tet8cube_face4_rank", my_rank, ".dat")
    # f = open(fname, "w")
    # for i = 1:length(mesh.vert_coords)
    #   println(f, mesh.vert_coords[i])
    # end
    # close(f)

    test_values = readdlm(fname)
    for i = 1:length(test_values)
      err = abs(test_values[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

  end # End context("Check control point manipulation on a geometry face")

  MPI.Barrier(comm)
  # Reset the coordinates and mesh to the original value
  for i = 1:mesh.numEl
    update_coords(mesh, i, orig_vert_coords[:,:,i])
  end
  commit_coords(mesh, sbp)

  context("--- Checking evaldXdControlPointProduct for 3D DG Mesh ---") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0.5, 0.5, 0.5] # No offset in the X & Y direction
    geom_faces = opts["BC2"]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

    fill!(map.work, 0.0)

    # Create seed vector
    # - Get original wall coordinates
    orig_wallCoords = FreeFormDeformation.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
    Xs_bar = ones(3, size(orig_wallCoords,2))
    evaldXdControlPointProduct(map, mesh, vec(Xs_bar))

    fname = "./testvalues/evaldXdControlPointProduct_tet8cube.dat"

    test_values = readdlm(fname)
    @fact length(map.work) --> length(test_values)
    for i = 1:length(map.work)
      err = abs(test_values[i] - map.work[i])
      @fact err --> less_than(1e-14) "problem at index $i"
    end

  end # End context("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")

end # End facts("--- Checking FFD on 3D serial DG Pumi meshes ---")
