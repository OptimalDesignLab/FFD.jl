# test components

facts("--- Checking Mapping object ---") do

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

  @fact map.ndim --> 3
  @fact map.nctl --> [3,3,3]
  @fact map.order --> [2,2,2]
  @fact map.numnodes --> [3,3,3]

  context("--- Checking Knot calculations ---") do

    calcKnot(map)
    for i = 1:map.ndim
      @fact map.edge_knot[i][1] --> 0.0
      @fact map.edge_knot[i][2] --> 0.0
      @fact map.edge_knot[i][3] --> roughly(0.5, atol=1e-15)
      @fact map.edge_knot[i][4] --> 1.0
      @fact map.edge_knot[i][5] --> 1.0
    end

  end  # End context("--- Checking Knot calculations ---")

end  # End facts("--- Checking Mapping object ---")

# Pumi Specific Tests
facts("--- Checking FFD Types and Functions For Serial DG Pumi Meshes ---") do

  # MPI Declarations
  MPI.Init()
  comm = MPI.COMM_WORLD
  comm_world = MPI.MPI_COMM_WORLD
  comm_self = MPI.COMM_SELF
  my_rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)

  opts = PdePumiInterface.get_defaults()
  # 2D mesh
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = true
  opts["operator_type"] = "SBPOmega"
  opts["dmg_name"] = "../src/mesh_files/2D_Airfoil.dmg"
  opts["smb_name"] = "../src/mesh_files/2D_Airfoil.smb"
  opts["numBC"] = 2

  # For 2DAirfoil
  opts["BC1"] = [8,11,14,17]
  opts["BC1_name"] = "FarField"
  opts["BC2"] = [5]
  opts["BC2_name"] = "Airfoil"

  opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
  opts["jac_type"] = 2

  # Create PumiMesh and SBP objects
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

  context("--- Checking Linear Mapping For Entire DG Mesh ---") do

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh)

    # Create knot vector
    calcKnot(map)

    # Create Bounding box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

    # Control points
    controlPoint(map, box)

    # Populate map.xi
    calcParametricMappingLinear(map, box, mesh)


    @fact map.ndim --> 2
    @fact map.full_geom --> true
    @fact map.nctl --> [4,4,2]
    @fact map.order --> [4,4,2]
    @fact map.edge_knot[1] --> [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @fact map.edge_knot[2] --> [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @fact map.edge_knot[3] --> [0.0,0.0,1.0,1.0]
    @fact size(map.aj) --> (3,4,3)
    @fact size(map.dl) --> (3,3)
    @fact size(map.dr) --> (3,3)
    @fact size(map.work) --> (4,4,2,12)
    @fact size(map.cp_xyz) --> (3,4,4,2)
    @fact size(map.xi) --> (3,3,1498)

    # Check Control Point Coordinates
    outname = string("./testvalues/control_points_full_linear_mapping.dat")
    orig_control_pts = readdlm(outname)
    for i = 1:length(map.cp_xyz)
      err = norm(map.cp_xyz[i] - orig_control_pts[i], 2)
      @fact err --> less_than(1e-14)
    end

    # Check map.xi
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @fact err --> less_than(1e-14)
    end

    # Bounding Box FactCheck
    @fact box.ndim --> 3
    @fact box.offset --> [0.0,0.0,0.5]
    @fact box.origin --> [-0.1,-0.1,-0.5]
    @fact box.unitVector --> [1.0 0.0 0.0
                              0.0 1.0 0.0
                              0.0 0.0 1.0]
    @fact box.geom_bounds --> [-0.1 -0.1 -0.0
                                0.2 0.1 0.0]
    @fact box.box_bounds --> [-0.1 -0.1 -0.5
                               0.2 0.1 0.5]
    @fact box.lengths[1] --> roughly(0.3, atol = 1e-14)
    @fact box.lengths[2] --> roughly(0.2, atol = 1e-14)
    @fact box.lengths[3] --> roughly(1.0, atol = 1e-14)

  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")

  # Free Form deformation parameters
  ndim = 2
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [4,4,2]

  # Create Mapping object
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, box)

  # Populate map.xi
  calcParametricMappingNonlinear(map, box, mesh)

  context("--- Checking Nonlinear Mapping For Entire DG Mesh ---") do

    # Check if the values are correct
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @fact err --> less_than(1e-14)
    end

  end # End context("--- Checking Nonlinear mapping For Entire DG Mesh ---")

  context("--- Checking Control Point Manipulation ---") do

    # Rigid body rotation
    theta = -20*pi/180  # Rotate wall coordinates by 10 degrees
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

    evalVolume(map, mesh)

    outname = string("./testvalues/volume_coords_full_DG_mesh_2D_airfoil.dat")
    test_vert_coords = readdlm(outname)
    @fact length(test_vert_coords) --> length(mesh.vert_coords)
    for i = 1:20:length(mesh.vert_coords)
      err = norm(test_vert_coords[i] - mesh.vert_coords[i], 2)
      @fact err --> less_than(1e-14)
    end

    for i = 1:mesh.numEl
      update_coords(mesh, i, mesh.vert_coords[:,:,i])
    end
    commit_coords(mesh, sbp)
    writeVisFiles(mesh, "translation_plus_rotation_DG")

  end


end # End facts("--- Checking FFD Types and Functions For Serial DG Pumi Meshes ---")

facts("--- Checking Functions Specific to CG Pumi Meshes in Serial ---") do

  # MPI Declarations
  comm = MPI.COMM_WORLD
  comm_world = MPI.MPI_COMM_WORLD
  comm_self = MPI.COMM_SELF
  my_rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)

  opts = PdePumiInterface.get_defaults()
  # 2D mesh
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = false
  opts["operator_type"] = "SBPGamma"
  opts["dmg_name"] = "../src/mesh_files/2D_Airfoil.dmg"
  opts["smb_name"] = "../src/mesh_files/2D_Airfoil.smb"
  opts["numBC"] = 2

  # For 2DAirfoil
  opts["BC1"] = [8,11,14,17]
  opts["BC1_name"] = "FarField"
  opts["BC2"] = [5]
  opts["BC2_name"] = "Airfoil"

  opts["coloring_distance"] = 0 # 0 For CG Mesh 2 for DG Mesh
  opts["jac_type"] = 2

  # Create PumiMesh and SBP objects
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

  context("--- Checking Linear Mapping For Entire CG Mesh ---") do

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh)

    # Create knot vector
    calcKnot(map) # No need to test this since doesnt depend on mesh type

    # Create Bounding box
    # Doesn't depend on mesh type if entire geometry is embedded in FFD box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

    # Control points
    controlPoint(map, box) # No Need to test this since it doesnt depend on mesh type

    # Populate map.xi
    calcParametricMappingLinear(map, box, mesh)

    # Check map.xi
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @fact err --> less_than(1e-14)
    end

  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")

  context("--- Checking Nonlinear Mapping For Entire CG Mesh ---") do

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh)

    # Create knot vector
    calcKnot(map) # No need to test this since doesnt depend on mesh type

    # Create Bounding box
    # Doesn't depend on mesh type if entire geometry is embedded in FFD box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

    # Control points
    controlPoint(map, box) # No Need to test this since it doesnt depend on mesh type

    # Populate map.xi
    calcParametricMappingNonlinear(map, box, mesh)

    # Check map.xi
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @fact err --> less_than(1e-14)
    end

  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")

end # End facts("--- Checking Functions Specific to CG Pumi Meshes in Serial ---")

MPI.Finalize()
#=
facts("--- Checking BoundingBox ---") do

  @fact box.ndim --> 3
  @fact box.origin --> [0.5, 0.5, 0.5]
  @fact box.unitVector --> [1. 0. 0.;0. 1. 0.;0. 0. 1.]
  @fact box.geom_bound --> [1. 1. 1.;3. 3. 3.]
  @fact box.offset --> [0.5,0.5,0.5]
  @fact box.box_bound --> [0.5 0.5 0.5;3.5 3.5 3.5]

end  # End facts("--- Checking BoundingBox ---")

facts("--- Checking Control Point Generation ---") do

  controlPoint(map, box)
  for i = 1:3
    @fact map.cp_xyz[1,1,1,i] --> roughly(0.5, atol = 1e-15)
    @fact map.cp_xyz[2,2,2,i] --> roughly(2.0, atol = 1e-15)
    @fact map.cp_xyz[3,3,3,i] --> roughly(3.5, atol = 1e-15)
  end
  @fact map.cp_xyz[1,2,3,1] --> roughly(0.5, atol = 1e-15)
  @fact map.cp_xyz[1,2,3,2] --> roughly(2.0, atol = 1e-15)
  @fact map.cp_xyz[1,2,3,3] --> roughly(3.5, atol = 1e-15)

end # End facts("--- Checking Contol Point Generation ---")

facts("--- Checking Linear Mapping ---") do

  calcParametricMappingLinear(map, box, nodes_xyz)
  for i = 1:3
    @fact map.xi[1,1,1,i] --> roughly(0.16666666666666666, atol = 1e-15)
    @fact map.xi[2,2,2,i] --> roughly(0.5, atol = 1e-15)
    @fact map.xi[3,3,3,i] --> roughly(0.8333333333333334, atol = 1e-15)
  end
  @fact map.xi[1,2,3,1] --> roughly(0.16666666666666666, atol = 1e-15)
  @fact map.xi[1,2,3,2] --> roughly(0.5, atol = 1e-15)
  @fact map.xi[1,2,3,3] --> roughly(0.8333333333333334, atol = 1e-15)

end # End facts("--- Checking Linear Mapping ---")

facts("--- Checking Nonlinear Mapping ---") do

  context("Checking Volume Derivative") do

    xi = [0.5,0.5,0.5]
    dX = zeros(AbstractFloat, 3)
    jderiv = [1,0,0]
    FreeFormDeformation.calcdXdxi(map, xi, jderiv, dX)
    @fact dX[1] --> roughly(3.0, atol = 1e-15)
    @fact dX[2] --> roughly(0.0, atol = 1e-15)
    @fact dX[3] --> roughly(0.0, atol = 1e-15)

  end  # End context("Checking Volume Derivative")

  context("Checking nonlinearMap") do

    X = ones(AbstractFloat, 3)
    pX = zeros(X)

    FreeFormDeformation.nonlinearMap(map, box, X, pX)
    for i = 1:3
      @fact pX[i] --> roughly(0.16666666666666666, atol = 1e-15)
    end

  end  # End context ("Check nonlinearMap")

  context("Checking calcParametricMappingNonlinear") do

    calcParametricMappingNonlinear(map, box, nodes_xyz)
    for i = 1:3
      @fact map.xi[1,1,1,i] --> roughly(0.16666666666666666, atol = 1e-15)
      @fact map.xi[2,2,2,i] --> roughly(0.5, atol = 1e-15)
      @fact map.xi[3,3,3,i] --> roughly(0.8333333333333334, atol = 1e-15)
    end
    @fact map.xi[1,2,3,1] --> roughly(0.16666666666666666, atol = 1e-15)
    @fact map.xi[1,2,3,2] --> roughly(0.5, atol = 1e-15)
    @fact map.xi[1,2,3,3] --> roughly(0.8333333333333334, atol = 1e-15)

  end  # End context ("Checking calcParametricMappingNonlinear")


end # End facts("--- Checking Linear Mapping ---")


facts("--- Checking FFD Volume Evaluation ---") do

  context("Checking single point evaluation") do

    xyz = zeros(AbstractFloat, map.ndim)
    xi = 0.5*ones(AbstractFloat, map.ndim)
    FreeFormDeformation.evalVolumePoint(map, xi, xyz)
    for i = 1:map.ndim
      @fact xyz[i] --> roughly(2.0, atol = 1e-15)
    end

  end  # End context("Checking single poitn evaluation")

  context("Checking multiple point evaluation") do
    Vol = zeros(nodes_xyz)
    FreeFormDeformation.evalVolume(map, Vol)
    for idim = 1:3
      for k = 1:map.numnodes[3]
        for j = 1:map.numnodes[2]
          for i = 1:map.numnodes[1]
            err = Vol[i,j,k,idim] - nodes_xyz[i,j,k,idim]
            @fact err --> roughly(0.0, atol = 1e-14)
          end
        end
      end
    end

  end # End context("Checking multiple point evaluation")

end  # facts("--- Checking FFD Volume Evaluation ---")
=#
