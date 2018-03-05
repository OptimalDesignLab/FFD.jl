# All 3D Pumi mesh tests go here

# MPI Declarations
opts = Dict{ASCIIString, Any}()
opts["order"] = 1
opts["dimensions"] = 3
opts["use_DG"] = true
opts["Tsbp"] = Float64
opts["Tmsh"] = Complex128
opts["operator_type"] = "SBPOmega"
opts["smb_name"] = "../src/mesh_files/tet8cube.smb"
opts["dmg_name"] = ".null"

opts["numBC"] = 3
opts["BC1"] = [ 0, 1, 2, 3]
opts["BC1_name"] = "ExpBC"
opts["BC2"] = [4]
opts["BC2_name"] = "noPenetrationBC"
opts["BC3"] = [10]
opts["BC3_name"] = "noPenetrationBC"

opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
opts["jac_type"] = 2
opts["jac_method"] = 2
opts["run_type"] = 5

Tsbp = opts["Tsbp"]
Tmsh = opts["Tmsh"]
shape_type = 2
sbp = getTetSBPOmega(degree=opts["order"], Tsbp=Tsbp)
ref_verts = sbp.vtx
face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
topo = ElementTopology{3}(face_verts)
sbpface = TetFace{Tsbp}(opts["order"], sbp.cub, ref_verts)
mesh = PumiMeshDG3(Tmsh, sbp, opts, sbpface, topo, dofpernode=dofpernode,
                     shape_type=shape_type)

#mesh = PumiMeshDG3{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"], sbp,
#                         opts, sbpface, topo; dofpernode=1,
#                         coloring_distance=opts["coloring_distance"],
#                         shape_type=shape_type)

orig_vert_coords = deepcopy(mesh.vert_coords)

facts("--- Checking FFD on 3D serial DG Pumi meshes ---") do

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
    commitToPumi(map, mesh, sbp, vertices, opts)

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
    update_coords(mesh, i, real(orig_vert_coords[:,:,i]))
  end
  commit_coords(mesh, sbp, opts)


  context("Check control point manipulation on a geometry face") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0.1, 0.1, 0.1] # No offset in the X & Y direction
#    geom_faces = opts["BC2"]
    bc_nums = [2]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, bc_nums)
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
    commitToPumi(map, mesh, sbp, vertices, opts)

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

  # Reset the coordinates and mesh to the original value
  for i = 1:mesh.numEl
    update_coords(mesh, i, real(orig_vert_coords[:,:,i]))
  end
  commit_coords(mesh, sbp, opts)

  context("--- Checking evaldXdControlPointProduct for 3D DG Mesh ---") do

    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0.5, 0.5, 0.5] # No offset in the X & Y direction
#    geom_faces = opts["BC2"]
    bc_nums = [2]

    # Create a mapping object using nonlinear mapping
    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, false, bc_nums)

    fill!(map.work, 0.0)
    pert = 1e-6

    # Create seed vector
    # - Get original wall coordinates
    orig_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, bc_nums)
    nwall_faces = FreeFormDeformation.getnWallFaces(mesh, bc_nums)
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

    # fname = "./testvalues/evaldXdControlPointProduct_tet8cube.dat"
    # f = open(fname, "w")
    # for i = 1:length(map.work)
    #   println(f, map.work[i])
    # end
    # close(f)

    # Check against finite difference
    cp_jacobian = zeros(length(orig_wallCoords), length(map.cp_xyz))
    for i = 1:length(map.cp_xyz)
      map.cp_xyz[i] += pert
      vertices = evalSurface(map, mesh)
      commitToPumi(map, mesh, sbp, vertices, opts)
      new_wallCoords = FreeFormDeformation.getUniqueWallCoordsArray(mesh, bc_nums)
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
