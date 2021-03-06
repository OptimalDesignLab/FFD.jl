# test components

@testset "--- Checking Generic Mapping object ---" begin

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

  @test ( map.ndim )== 3
  @test ( map.nctl )== [3,3,3]
  @test ( map.order )== [2,2,2]
  @test ( map.numnodes )== [3,3,3]

  @testset "--- Checking Knot calculations ---" begin

    calcKnot(map)
    for i = 1:map.ndim
      @test ( map.edge_knot[i][1] )== 0.0
      @test ( map.edge_knot[i][2] )== 0.0
      @test isapprox( map.edge_knot[i][3], 0.5) atol=1e-15
      @test ( map.edge_knot[i][4] )== 1.0
      @test ( map.edge_knot[i][5] )== 1.0
    end
  end  # End context("--- Checking Knot calculations ---")

end  # End facts("--- Checking Mapping object ---")

# MPI Declarations
comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

# Pumi Specific Tests
@testset "--- Checking FFD Types and Functions For Full Serial DG Pumi Meshes ---" begin

  # 2D mesh
  opts = Dict{ASCIIString, Any}()
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = true
  opts["Tsbp"] = Float64
  opts["Tmsh"] = Complex128
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
  opts["jac_method"] = 2
  opts["run_type"] = 5

  # Create PumiMesh and SBP objects
  dofpernode = 1
  ref_verts = [-1. 1 -1; -1 -1 1]
  Tsbp = opts["Tsbp"]
  Tmsh = opts["Tmsh"]
  sbp = getTriSBPOmega(degree=opts["order"], Tsbp=opts["Tsbp"])
  sbpface = TriFace{Tsbp}(opts["order"], sbp.cub, ref_verts.')
  topo = 0
  shape_type = 2
  mesh = PumiMeshDG2(Tmsh, sbp, opts, sbpface, dofpernode=dofpernode,
                     shape_type=shape_type)
#  mesh = PumiMeshDG2{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"],
#                           sbp, opts, sbpface; dofpernode=dofpernode,
#                           coloring_distance=opts["coloring_distance"],
#                           shape_type=shape_type)

  @testset "--- Checking Linear Mapping For DG Mesh ---" begin

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, bc_nums=[1, 2])

    # Create knot vector
    calcKnot(map)

    # Create Bounding box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)
#=
    # Control points
    controlPoint(map, box)

    # Populate map.xi
    calcParametricMappingLinear(map, box, mesh, )


    @test ( map.ndim )== 2
    @test ( map.nctl )== [4,4,2]
    @test ( map.order )== [4,4,2]
    @test ( map.edge_knot[1] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[2] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[3] )== [0.0,0.0,1.0,1.0]
    @test ( size(map.aj) )== (3,4,3)
    @test ( size(map.dl) )== (3,3)
    @test ( size(map.dr) )== (3,3)
    @test ( size(map.work) )== (12,4,4,2)
    @test ( size(map.cp_xyz) )== (3,4,4,2)
    @test ( size(map.xi) )== (3,3,1498)
=#
#=
    # Check Control Point Coordinates
    outname = string("./testvalues/control_points_full_linear_mapping.dat")
    orig_control_pts = readdlm(outname)
    for i = 1:length(map.cp_xyz)
      err = norm(map.cp_xyz[i] - orig_control_pts[i], 2)
      @test  err  < 1e-14
    end
=#
#=
    # Check map.xi
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @test  err  < 1e-14
    end
=#
    # Bounding Box Base.Test
    @test ( box.ndim )== 3
    @test ( box.offset )== [0.0,0.0,0.5]
    @test ( box.origin )== [-0.1,-0.1,-0.5]
    @test ( box.unitVector )== [1.0 0.0 0.0
                              0.0 1.0 0.0
                              0.0 0.0 1.0]
    @test ( box.geom_bounds )== [-0.1 -0.1 -0.0
                                0.2 0.1 0.0]
    @test ( box.box_bounds )== [-0.1 -0.1 -0.5
                               0.2 0.1 0.5]
    @test isapprox( box.lengths[1], 0.3) atol= 1e-14
    @test isapprox( box.lengths[2], 0.2) atol= 1e-14
    @test isapprox( box.lengths[3], 1.0) atol= 1e-14
  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")
end  # end facts block
#=
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

  @testset "--- Checking Nonlinear Mapping For Entire DG Mesh ---" begin
#=
    # Check if the values are correct
    outname = string("./testvalues/xi_full_linear_mapping_DG.dat")
    orig_xi_vals = readdlm(outname)
    for i = 1:10:length(map.xi)
      err = norm(orig_xi_vals[i] - map.xi[i], 2)
      @test  err  < 1e-14
    end
=#
  end # End context("--- Checking Nonlinear mapping For Entire DG Mesh ---")

  @testset "--- Checking Control Point Manipulation ---" begin
#=
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

    vertices = evalVolume(map, mesh)
    commitToPumi(map, mesh, sbp, vertices, opts)

    outname = string("./testvalues/volume_coords_full_DG_mesh_2D_airfoil.dat")
    test_vert_coords = readdlm(outname)
    @test ( length(test_vert_coords) )== length(mesh.vert_coords)
    for i = 1:20:length(mesh.vert_coords)
      err = norm(test_vert_coords[i] - mesh.vert_coords[i], 2)
      @test  err  < 1e-14
    end

    writeVisFiles(mesh, "translation_plus_rotation_DG")
=#
  end


end # End facts("--- Checking FFD Types and Functions For Serial DG Pumi Meshes ---")

@testset "--- Checking Specific Geometry Faces in Pumi DG Mesh Embedded in FFD ---" begin

  comm = MPI.COMM_WORLD
  comm_world = MPI.MPI_COMM_WORLD
  comm_self = MPI.COMM_SELF
  my_rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)

  opts = Dict{ASCIIString, Any}()
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = true
  opts["Tsbp"] = Float64
  opts["Tmsh"] = Complex128
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
  opts["jac_method"] = 2
  opts["run_type"] = 5

  # Create PumiMesh and SBP objects
  dofpernode = 1
  ref_verts = [-1. 1 -1; -1 -1 1]
  Tsbp = opts["Tsbp"]
  Tmsh = opts["Tmsh"]
  sbp = getTriSBPOmega(degree=opts["order"], Tsbp=opts["Tsbp"])
  sbpface = TriFace{Tsbp}(opts["order"], sbp.cub, ref_verts.')
  topo = 0
  shape_type = 2
  mesh = PumiMeshDG2(Tmsh, sbp, opts, sbpface, dofpernode=dofpernode,
                     shape_type=shape_type)

#  mesh = PumiMeshDG2{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"],
#                           sbp, opts, sbpface; dofpernode=dofpernode,
#                           coloring_distance=opts["coloring_distance"],
#                           shape_type=shape_type)

  @testset "--- Checking Linear Mapping for 2D DG Mesh ---" begin

    # geometry faces to be embedded in FFD Box
#    geom_faces = opts["BC2"]
    bc_nums =[2]

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, bc_nums=bc_nums)

    # Create knot vector
    calcKnot(map)

    # Create Bounding box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

    # Control points
    controlPoint(map, box)
#=
    calcParametricMappingLinear(map, box, mesh, bc_nums)

    @test ( map.ndim )== 2
    @test ( map.full_geom )== false
    @test ( map.nctl )== [4,4,2]
    @test ( map.order )== [4,4,2]
    @test ( map.edge_knot[1] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[2] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[3] )== [0.0,0.0,1.0,1.0]
    @test ( size(map.aj) )== (3,4,3)
    @test ( size(map.dl) )== (3,3)
    @test ( size(map.dr) )== (3,3)
    @test ( size(map.work) )== (12,4,4,2)
    @test ( size(map.cp_xyz) )== (3,4,4,2)
    @test ( size(map.xi) )== (1,)
    @test ( size(map.xi[1]) )== (3,2,102) # Essentially tests defineMapXi
=#
#=
    outname = string("./testvalues/xi_values_2D_airfoil_face5.dat")
    ctr = 1
    test_xi_values = readdlm(outname)
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:end_index
      bndry_facenums =sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(map.xi[itr][1,j,i] - test_xi_values[ctr],2)
          err2 = norm(map.xi[itr][2,j,i] - test_xi_values[ctr+1],2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)
=#

  end # End context("--- Checking Linear Mapping for DG Mesh ---")

  # geometry faces to be embedded in FFD Box
#  geom_faces = opts["BC2"]
  bc_nums =[2]

  # Free Form deformation parameters
  ndim = 2
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [4,4,2]

  # Create Mapping object
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, bc_nums=bc_nums)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, box)

  @testset "--- Checking Nonlinear Mapping for DG Mesh ---" begin
#=
    # Populate map.xi
    calcParametricMappingNonlinear(map, box, mesh, geom_faces)
    @test ( size(map.xi) )== (1,)
    @test ( size(map.xi[1]) )== (3,2,102) # Essentially tests defineMapXi

    outname = string("./testvalues/xi_values_2D_airfoil_face5.dat")
    ctr = 1
    test_xi_values = readdlm(outname)
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:end_index
      bndry_facenums =sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(map.xi[itr][1,j,i] - test_xi_values[ctr],2)
          err2 = norm(map.xi[itr][2,j,i] - test_xi_values[ctr+1],2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)
=#
  end # End context("--- Checking Nonlinear Mapping for DG Mesh ---")

  @testset "--- Checking Surface evaluation for 2D DG Mesh ---" begin

#=
#    commitToPumi(map, mesh, sbp, vertices, opts)
#    writeVisFiles(mesh, "mesh_warped")

    outname = string("./testvalues/modified_coordinates_2D_airfoil_face5.dat")
    test_surface_coords = readdlm(outname)
    ctr = 1
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums =sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(mesh.vert_coords[1,vtx_arr[j],bndry_i.element] -
                 test_surface_coords[ctr], 2)
          err2 = norm(mesh.vert_coords[2,vtx_arr[j],bndry_i.element] -
                 test_surface_coords[ctr+1], 2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
=#
  end # End context("--- Checking Surface evaluation for 2D DG Mesh ---")
#=
  @testset "--- Checking evaldXdControlPointProduct for 2D DG Mesh ---" begin

    fill!(map.work, 0.0)
    pert = 1e-6

    # Create seed vector
    # - Get original wall coordinates
    orig_wallCoords = evalSurface(map, mesh)[1]  # only do the first geometric entity
#    orig_wallCoords = FFD.getUniqueWallCoordsArray(mesh, geom_faces)
    nwall_faces = FFD.getnWallFaces(mesh, geom_faces)
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
      commitToPumi(map, mesh, sbp, vertices, opts)
      new_wallCoords = FFD.getUniqueWallCoordsArray(mesh, geom_faces)
      cp_jacobian[:,i] = (vec(new_wallCoords) - vec(orig_wallCoords))/pert
      map.cp_xyz[i] -= pert
    end # End for i = 1:length(map.cp_xyz)
    prod_val = transpose(cp_jacobian)*vec(Xs_bar)

    error = vec(cp_xyz_bar) - prod_val
    for i = 1:length(error)
      @test isapprox( error[i], 0.0) atol=1e-9
    end

  end # End context("--- Checking evaldXdControlPointProduct for 2D DG Mesh ---")
=#
  @testset "---- Checking Transposed Jacobian-vector Product -----" begin
#=
    delta_B = zeros(eltype(map.work), 3, map.nctl[1], map.nctl[2], map.nctl[3])
    Xs_orig = evalSurface(map, mesh)
    println("typeof(Xs_orig) = ", typeof(Xs_orig))
    println("size(Xs_orig) = ", size(Xs_orig))
    println("size(Xs_orig[1]) = ", size(Xs_orig[1]))
      
    delta_S = Array{Array{Complex128, 3}}(length(map.geom_faces))

    h = 1e-6
    
    for itr=1:length(map.geom_faces)
      println("itr = ", itr)
      delta_S[itr] = zeros(Complex128, size(Xs_orig[itr]))

      # compute finite differences delta_s^T jac * delta_ba
      # and compare against the AD version
      # here we take all combinations of delta_S and delta_x = 1 at a single
      # entry
      for i=4:length(map.cp_xyz)
        # in 2D, skip z components
        if (i % 3) == 0
          continue
        end
        println("i = ", i)
        map.cp_xyz[i] += h
        verts_new = evalSurface(map, mesh)
        map.cp_xyz[i] -= h

        for j=1:length(verts_new[itr])
#          if (j % 3) == 0
#            continue
#          end
          println("j = ", j)
          fill!(delta_B, 0.0)
          delta_S[itr][j] = 1
          FFD.evaldXdControlPointTransposeProduct(map, mesh, delta_S, delta_B)

          delta_xj = (verts_new[itr][j] - Xs_orig[itr][j])/h
          println("delta_xj = ", delta_xj)
          println("delta_B = ", delta_B[i])

          @test isapprox( delta_B[i], delta_xj) atol=1e-6
          delta_S[itr][j] = 0
        end
      end
    end
=#

      end  # end context

end # End facts("--- Checking Specific Geometry Faces in Pumi DG Mesh Embedded in FFD ---")

#=
@testset "--- Checking Functions Specific to CG Pumi Meshes in Serial ---" begin

  # MPI Declarations
  comm = MPI.COMM_WORLD
  comm_world = MPI.MPI_COMM_WORLD
  comm_self = MPI.COMM_SELF
  my_rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)

  opts = Dict{ASCIIString, Any}()
  opts["order"] = 1
  opts["dimensions"] = 2
  opts["use_DG"] = true
  opts["Tsbp"] = Float64
  opts["Tmsh"] = Complex128
  opts["operator_type"] = "SBPGamma"
  opts["dmg_name"] = "../src/mesh_files/2D_Airfoil.dmg"
  opts["smb_name"] = "../src/mesh_files/2D_Airfoil.smb"
  opts["numBC"] = 2

  # For 2DAirfoil
  opts["BC1"] = [8,11,14,17]
  opts["BC1_name"] = "FarField"
  opts["BC2"] = [5]
  opts["BC2_name"] = "Airfoil"

  opts["coloring_distance"] = 0 # For CG Mesh 2 for DG Mesh
  opts["jac_type"] = 2
  opts["jac_method"] = 2
  opts["run_type"] = 5
  opts["use_edge_res"] = false
  opts["write_edge_vertnums"] = false
  opts["write_face_vertnums"] = false

  # Create PumiMesh and SBP objects
  dofpernode = 1
  ref_verts = [-1. 1 -1; -1 -1 1]
  Tsbp = opts["Tsbp"]
  Tmsh = opts["Tmsh"]

  shape_type = 1
  sbp = getTriSBPGamma(degree=opts["order"], Tsbp=Tsbp)
  sbpface = TriFace{Float64}(opts["order"], sbp.cub, sbp.vtx)
  mesh = PumiMesh2{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"], sbp, opts, sbpface;
                             dofpernode=dofpernode,
                             coloring_distance=opts["coloring_distance"],
                             shape_type=shape_type)

  @testset "--- Checking Linear Mapping For Entire CG Mesh ---" begin

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
      @test  err  < 1e-14
    end

  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")

  @testset "--- Checking Nonlinear Mapping For Entire CG Mesh ---" begin

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
      @test  err  < 1e-14
    end

  end # End context("--- Checking Linear mapping For Entire DG Mesh ---")

end # End facts("--- Checking Functions Specific to CG Pumi Meshes in Serial ---")

@testset "--- Checking Specific Geometry Faces in Pumi CG Mesh Embedded in FFD ---" begin

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
  opts["jac_method"] = 2
  opts["run_type"] = 5

  # Create PumiMesh and SBP objects
  sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

  @testset "--- Checking Linear Mapping for 2D CG Mesh ---" begin

    # geometry faces to be embedded in FFD Box
    geom_faces = opts["BC2"]

    # Free Form deformation parameters
    ndim = 2
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]

    # Create Mapping object
    map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, geom_faces=geom_faces)

    # Create knot vector
    calcKnot(map)

    # Create Bounding box
    offset = [0., 0., 0.5] # No offset in the X & Y direction
    box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

    # Control points
    controlPoint(map, box)

    calcParametricMappingLinear(map, box, mesh, geom_faces)

    @test ( map.ndim )== 2
    @test ( map.full_geom )== false
    @test ( map.nctl )== [4,4,2]
    @test ( map.order )== [4,4,2]
    @test ( map.edge_knot[1] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[2] )== [0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0]
    @test ( map.edge_knot[3] )== [0.0,0.0,1.0,1.0]
    @test ( size(map.aj) )== (3,4,3)
    @test ( size(map.dl) )== (3,3)
    @test ( size(map.dr) )== (3,3)
    @test ( size(map.work) )== (12,4,4,2)
    @test ( size(map.cp_xyz) )== (3,4,4,2)

    outname = string("./testvalues/xi_values_2D_airfoil_face5.dat")
    test_xi_values = readdlm(outname)
    ctr = 1
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:end_index
      bndry_facenums =sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(map.xi[itr][1,j,i] - test_xi_values[ctr],2)
          err2 = norm(map.xi[itr][2,j,i] - test_xi_values[ctr+1],2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)

  end # End context("--- Checking Linear Mapping for DG Mesh ---")

  # geometry faces to be embedded in FFD Box
  geom_faces = opts["BC2"]

  # Free Form deformation parameters
  ndim = 2
  order = [4,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [4,4,2]

  # Create Mapping object
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, geom_faces=geom_faces)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  offset = [0., 0., 0.5] # No offset in the X & Y direction
  box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, box)

  @testset "--- Checking Nonlinear Mapping for CG Mesh ---" begin

    # Populate map.xi
    calcParametricMappingNonlinear(map, box, mesh, geom_faces)

    # Check for errors in map.xi
    outname = string("./testvalues/xi_values_2D_airfoil_face5.dat")
    test_xi_values = readdlm(outname)
    ctr = 1
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:end_index
      bndry_facenums =sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(map.xi[itr][1,j,i] - test_xi_values[ctr],2)
          err2 = norm(map.xi[itr][2,j,i] - test_xi_values[ctr+1],2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)

  end # End context("--- Checking Nonlinear Mapping for DG Mesh ---")

  @testset "--- Checking Surface evaluation for 2D CG Mesh ---" begin

    evalSurface(map, mesh, sbp)
    outname = string("./testvalues/modified_coordinates_2D_airfoil_face5.dat")
    #=
    test_surface_coords = readdlm(outname)
    ctr = 1
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums =sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          err1 = norm(mesh.coords[1,vtx_arr[j],bndry_i.element] -
                 test_surface_coords[ctr], 2)
          err2 = norm(mesh.coords[2,vtx_arr[j],bndry_i.element] -
                 test_surface_coords[ctr+1], 2)
          @test  err1  < 1e-14
          @test  err2  < 1e-14
          ctr += 2
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
    =#

  end # End context("--- Checking Nonlinear Mapping for DG Mesh ---")
end # End facts("--- Checking Specific Geometry Faces in Pumi CG Mesh Embedded in FFD ---")
=#
=#
