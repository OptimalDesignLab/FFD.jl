# All 3D Pumi mesh tests go here

# MPI Declarations
if !MPI.Initialized()
  MPI.Init()
end

resize!(ARGS, 1)
ARGS[1] = "./input_vals_3d.jl"

opts = PDESolver.read_input(ARGS[1])
sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)
orig_vert_coords = deepcopy(mesh.vert_coords)

facts("--- Checking FFD on 3D serial DG Pumi meshes ---") do
#=
  context("Check control point manipulation with nonlinear mapping on full mesh") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0., 0., 0.] # No offset in the X & Y direction

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, true)
    writeControlPointsVTS(map)

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

    writeVisFiles(mesh, "FFD_full_body_3D")

    fname = "./testvalues/translation_plus_rotation_DG_3D_tet8cube.dat"
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

  # Reset the coordinates and mesh to the original value
  for i = 1:mesh.numEl
    update_coords(mesh, i, orig_vert_coords[:,:,i])
  end
  commit_coords(mesh, sbp)


  context("Check control point manipulation on a geometry face") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0.1, 0.1, 0.1] # No offset in the X & Y direction
    geom_faces = opts["BC2"]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, geom_faces)
    writeControlPointsVTS(map)

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

    writeVisFiles(mesh, "FFD_perturbed_face_BC2")

    fname = "./testvalues/translation_plus_rotation_DG_3D_tet8cube_face4.dat"
    # f = open(fname, "w")
    # for i = 1:length(mesh.vert_coords)
    #  println(f, mesh.vert_coords[i])
    # end
    # close(f)

    test_values = readdlm(fname)
    for i = 1:length(test_values)
      err = abs(test_values[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

  end # End context("Check control point manipulation on a geometry face")
=#
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
    pert = 1e-6

    # Create seed vector
    # - Get original wall coordinates
    orig_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, geom_faces)
    nwall_faces = FreeFormDeformation.getnWallFaces(mesh, geom_faces)
    # Xs_bar = randn(3, size(orig_wallCoords,2))
    Xs_bar = ones(3, size(orig_wallCoords,2))
    cp_xyz_bar = zeros(map.cp_xyz)
    evaldXdControlPointProduct(map, mesh, vec(Xs_bar))
    for i = 1:size(map.work, 4)
      for j = 1:size(map.work, 3)
        for k = 1:size(map.work, 2)
            cp_xyz_bar[1:3,k,j,i] = map.work[1:3, k,j,i]
        end
      end
    end

    fname = "./testvalues/evaldXdControlPointProduct_tet8cube.dat"
    f = open(fname, "w")
    for i = 1:length(map.work)
      println(f, map.work[i])
    end
    close(f)

    # Check against finite difference
    cp_jacobian = zeros(length(orig_wallCoords), length(map.cp_xyz))
    for i = 1:length(map.cp_xyz)
      map.cp_xyz[i] += pert
      vertices = evalSurface(map, mesh)
      commitToPumi(map, mesh, sbp, vertices)
      new_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, geom_faces)
      cp_jacobian[:,i] = (vec(new_wallCoords) - vec(orig_wallCoords))/pert
      map.cp_xyz[i] -= pert
    end # End for i = 1:length(map.cp_xyz)
    prod_val = transpose(cp_jacobian)*vec(Xs_bar)

    error = vec(cp_xyz_bar) - prod_val
    for i = 1:length(error)
      @fact error[i] --> roughly(0.0, atol=1e-8) "Error problem at i = $i"
    end

  end # End context("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")

end # End facts("--- Checking FFD on 3D serial DG Pumi meshes ---")
