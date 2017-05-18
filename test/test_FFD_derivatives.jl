# Test derivatives w.r.t control points
facts("---Checking contractWithdGdB ---") do

  opts = PdePumiInterface.get_defaults()
  # 2D mesh
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = true
  opts["operator_type"] = "SBPOmega"
  opts["dmg_name"] = ".null"
  opts["smb_name"] = "../src/mesh_files/abc.smb"
  opts["numBC"] = 2

  opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
  opts["jac_type"] = 2
  opts["jac_method"] = 2
  opts["run_type"] = 5

  opts["numBC"] = 4
  opts["BC1"] = [0]
  opts["BC1_name"] = "FreeStreamBC"
  opts["BC2"] = [1]
  opts["BC2_name"] = "FreeStreamBC"
  opts["BC3"] = [2]
  opts["BC3_name"] = "FreeStreamBC"
  opts["BC4"] = [3]
  opts["BC4_name"] = "noPenetrationBC"

  # Create PumiMesh and SBP objects
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)
  # geometry faces to be embedded in FFD Box
  geom_faces = opts["BC4"]

  # Free Form deformation parameters
  ndim = mesh.dim
  order = [2,2,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,2,2]

  # Create Mapping object
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, geom_faces=geom_faces)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, box)
  writeControlPointsVTS(map)

  calcParametricMappingNonlinear(map, box, mesh, geom_faces)

  pert = 1e-6 # Finite difference perturbation for tests
  xi = map.xi[1][:,2,1]
  dJdG = rand(Float64, 3)
  dJdG[3] = 0.0

  FreeFormDeformation.contractWithdGdB(map, xi, dJdG)
  cp_xyz_bar = zeros(map.cp_xyz) # Extract the values out of map.work
  for i = 1:size(map.work,4)
    for j = 1:size(map.work,3)
      for k = 1:size(map.work,2)
        cp_xyz_bar[1:3,k,j,i] = map.work[1:3,k,j,i]
      end
    end
  end

  # Check against finite difference
  # - Get original wall coordinates
  orig_wallCoords = getUniqueWallCoordsArray(mesh, geom_faces, false)
  multiplying_vec = zeros(length(orig_wallCoords))
  multiplying_vec[4:5] = dJdG[1:2]
  cp_jacobian = zeros(length(orig_wallCoords), length(map.cp_xyz))
  for i = 1:length(map.cp_xyz)
    map.cp_xyz[i] += pert
    vertices = evalSurface(map, mesh)
    commitToPumi(map, mesh, sbp, vertices)
    new_wallCoords = getUniqueWallCoordsArray(mesh, geom_faces, false)
    cp_jacobian[:,i] = (vec(new_wallCoords) - vec(orig_wallCoords))/pert
    map.cp_xyz[i] -= pert
  end # End for i = 1:length(map.cp_xyz)

  prod_val = transpose(cp_jacobian)*multiplying_vec

  # Compute the error between reverse mode and finite difference
  error = vec(cp_xyz_bar) - prod_val
  for i = 1:length(error)
    @fact error[i] --> roughly(0.0, atol=1e-10)
  end

end # End facts("---Checking contractWithdGdB ---")

#=
facts("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---") do

  fill!(map.work, 0.0)

  # Create seed vector
  # - Get original wall coordinates
  orig_wallCoords = getUniqueWallCoordsArray(mesh, geom_faces, false)
  println("orig_wallCoords = \n$(orig_wallCoords)")
  nwall_faces = getnWallFaces(mesh, geom_faces)
  Xs_bar = randn(3, size(orig_wallCoords,2))
  Xs_bar[3,:] = 0.0 # To accurately simulate a 2D mesh
  cp_xyz_bar = zeros(map.cp_xyz)
  evaldXdControlPointProduct(map, mesh, vec(Xs_bar))
  for i = 1:size(map.work, 4)
    for j = 1:size(map.work, 3)
      for k = 1:size(map.work, 2)
          cp_xyz_bar[1:3,k,j,i] = map.work[1:3, k,j,i]
      end
    end
  end


  # Check against finite difference
  cp_jacobian = zeros(length(orig_wallCoords), length(map.cp_xyz))
  for i = 1:length(map.cp_xyz)
    map.cp_xyz[i] += pert
    vertices = evalSurface(map, mesh)
    commitToPumi(map, mesh, sbp, vertices)
    new_wallCoords = getUniqueWallCoordsArray(mesh, geom_faces, false)
    cp_jacobian[:,i] = (vec(new_wallCoords) - vec(orig_wallCoords))/pert
    map.cp_xyz[i] -= pert
  end # End for i = 1:length(map.cp_xyz)
  prod_val = transpose(cp_jacobian)*vec(Xs_bar)
  println("prod_val = \n$prod_val")
  println("vec(cp_xyz_bar) = \n$(vec(cp_xyz_bar))")
  # println("vec(Xs_bar) = \n$(vec(Xs_bar))")

  error = vec(cp_xyz_bar) - prod_val
  for i = 1:length(error)
    # println("error[$i] = $(error[i])")
    @fact error[i] --> roughly(0.0, atol=1e-9) "Error problem at i = $i"
  end

end # End facts("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")
=#
