# Test FFD function that make parallel MPI calls


# MPI Declarations
comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

# opts = PdePumiInterface.get_defaults()
# 2D mesh
opts = Dict{ASCIIString, Any}()
opts["order"] = 1
opts["dimensions"] = 2
opts["use_DG"] = true
opts["operator_type"] = "SBPOmega"
opts["dmg_name"] = ".null"
opts["smb_name"] = "../src/mesh_files/gvortex1np2.smb"
opts["Tsbp"] = Float64
opts["Tmsh"] = Complex128

# For 2D Isentropic Vortex
opts["numBC"] = 4
opts["BC1"] = [0]
opts["BC1_name"] = "FreeStreamBC"
opts["BC2"] = [1]
opts["BC2_name"] = "FreeStreamBC"
opts["BC3"] = [2]
opts["BC3_name"] = "FreeStreamBC"
opts["BC4"] = [3]
opts["BC4_name"] = "noPenetrationBC"

opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
opts["jac_type"] = 1
opts["run_type"] = 1
opts["jac_method"] = 2
opts["use_jac_precond"] = true

# Create PumiMesh and SBP objects
Tsbp = opts["Tsbp"]
Tmsh = opts["Tmsh"]
sbp = getTriSBPOmega(degree=opts["order"], Tsbp=Tsbp)
shape_type = 2
ref_verts = [-1. 1 -1; -1 -1 1]
sbpface = TriFace{Tsbp}(opts["order"], sbp.cub, ref_verts.')
topo = 0
mesh = PumiMeshDG2{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"],
                         sbp, opts, sbpface; dofpernode=1,
                         coloring_distance=opts["coloring_distance"],
                         shape_type=shape_type)
# sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

orig_vert_coords = deepcopy(mesh.vert_coords)

facts("--- Checking control point manipulation for 2D parallel mesh") do

  context("Check for the entire mesh") do
    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0., 0., 0.5] # No offset in the X & Y direction

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, true)

    # Control Point Manipulation
    #  - Rigid body rotation
    theta = -90*pi/180  # Rotate wall coordinates by 90 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
    # - Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # - Rigid body translation
    map.cp_xyz[1,:,:,:] += 0.2
    map.cp_xyz[2,:,:,:] += 0.3

    vertices = evalVolume(map, mesh)
    commitToPumi(map, mesh, sbp, vertices, opts)

    filename = string("translation_plus_rotation_DG_parallel_fullmesh_", mesh.myrank, ".dat")
    # f = open(string("./testvalues/", filename), "w")
    # for i = 1:length(vertices)
    #   println(f, vertices[i])
    # end
    # close(f)

    test_vtx_coords = readdlm(string("./testvalues/", filename))
    for i = 1:length(test_vtx_coords)
      err = abs.(test_vtx_coords[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

    # writeVisFiles(mesh, "translation_plus_rotation_DG_parallel_fullmesh")

  end # End context("Check for the entire mesh")

  # Reset the coordinates and mesh to the original value
  for i = 1:mesh.numEl
    update_coords(mesh, i, real(orig_vert_coords[:,:,i]))
  end
  commit_coords(mesh, sbp, opts)

  context("Check for a geometric surface") do

    # Free Form deformation parameters
    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    Tmsh = Float64
    geom_faces = opts["BC4"]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

    # Control Point Manipulation
    #  - Rigid body rotation
    theta = -5*pi/180  # Rotate wall coordinates by 5 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
    # - Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # - Rigid body translation
    map.cp_xyz[1,:,:,:] += 0.05

    vertices = evalSurface(map, mesh)
    commitToPumi(map, mesh, sbp, vertices, opts)

    filename = string("translation_plus_rotation_DG_parallel_surfaceBC4_", mesh.myrank, ".dat")
    # f = open(string("./testvalues/", filename), "w")
    # for i = 1:length(mesh.vert_coords)
    #  println(f, mesh.vert_coords[i])
    # end
    # close(f)

    test_vtx_coords = readdlm(string("./testvalues/", filename))
    for i = 1:length(test_vtx_coords)
      err = abs.(test_vtx_coords[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

    writeVisFiles(mesh, "translation_plus_rotation_DG_parallel_surface")

  end # End context("Check for a geometric surface")

end # End facts("--- Checking control point manipulation for 2D parallel mesh")


# Reset the coordinates and mesh to the original value
for i = 1:mesh.numEl
  update_coords(mesh, i, real(orig_vert_coords[:,:,i]))
end
commit_coords(mesh, sbp, opts)
MPI.Barrier(comm)

facts("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---") do

  # Free Form deformation parameters
  ndim = mesh.dim
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [4,4,2]
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  Tmsh = Float64
  geom_faces = opts["BC4"]

  # Create a mapping object using nonlinear mapping
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)

  fill!(map.work, 0.0)

  # Create seed vector
  # - Get original wall coordinates
  orig_wallCoords = FFD.getGlobalUniqueWallCorrdsArray(mesh, geom_faces)
  Xs_bar = ones(3, size(orig_wallCoords,2))
  Xs_bar[3,:] = 0.0 # To accurately simulate a 2D mesh
  cp_xyz_bar = zeros(map.cp_xyz)
  evaldXdControlPointProduct(map, mesh, vec(Xs_bar))


  filename = "./testvalues/evaldXdControlPointProduct.dat"
  test_values = readdlm(filename)
  for i = 1:length(test_values)
    err = test_values[i] - map.work[i]
    @fact err --> less_than(1e-14)
  end


end # End
