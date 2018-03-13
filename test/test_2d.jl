# tests the Pumi interface to FFD

push!(LOAD_PATH, "../src/")

using FFD
using PumiInterface
using PdePumiInterface
using SummationByParts
using ODLCommonTools
using FactCheck


function getTestMesh{Tmsh}(dim::Integer, ::Type{Tmsh})

  degree = 1
  shape_type = 2


  if dim == 2
    sbp = getTriSBPOmega(degree=degree)
  sbpface = TriFace{Float64}(degree, sbp.cub, sbp.vtx)
  elseif dim == 3
    sbp = getTetSBPOmega(degree=degree)
    sbpface = TetFace{Float64}(degree, sbp.cub, sbp.vtx)
  end


  opts = PdePumiInterface.get_defaults()
  opts["smb_name"] = "meshes/square5.smb"
  opts["dmg_name"] = ".null"
  opts["order"] = 1
  opts["coloring_distance"] = 2
  opts["numBC"] = 2
  opts["BC1"] = [0]
  opts["BC2"] = [1,2,3]
  opts["BC1_name"] = "fakeBC"
  opts["BC2_name"] = "fake2BC"

  if dim == 3
    opts["dimensions"] = 3
    opts["BC2"] = [1,2,3,4,5]
    opts["smb_name"] = "meshes/tet4.smb"
    face_verts = SummationByParts.SymCubatures.getfacevertexindices(sbp.cub)
    topo = ElementTopology{3}(face_verts)
  end

  if dim == 2
    mesh = PumiMeshDG2(Tmsh, sbp, opts, sbpface, dofpernode=1,
                      shape_type=shape_type)
  else
    mesh = PumiMeshDG3(Tmsh, sbp, opts, sbpface, topo, dofpernode=1,
                      shape_type=shape_type)
  end


  return mesh, sbp, opts
end


function test_surface(map, mesh)


  facts("----- Testing FFD Surface Evaluation -----") do
    vertices_orig = zeros(Complex128, mesh.dim, map.numFacePts)
    evalSurface(map, vertices_orig)

    cp_xyz = getControlPoints(map)
    # Rigid body rotation
    theta = -20*pi/180  # Rotate wall coordinates by 10 degrees
    phi = -10*pi/180  # rotation about x axis (3D only)
    if mesh.dim == 2
      rotMat = [cos(theta) -sin(theta)
                sin(theta) cos(theta)] # Rotation matrix

      for i=1:map.nctl[2]
        for j=1:map.nctl[1]
          cp_xyz[:, j, i] = rotMat*cp_xyz[:, j, i]
        end
      end
    else
      rotMat1 = [cos(theta) -sin(theta) 0
                 sin(theta) cos(theta)  0
                 0          0           1] # Rotation matrix
      rotMat2 = [1 0         0        ;
                 0 cos(phi) -sin(phi) ;
                 0 sin(phi)  cos(phi) ]

      rotMat = rotMat1*rotMat2
      # Rotate the control points
      for k = 1:map.nctl[3]
        for j = 1:map.nctl[2]
          for i = 1:map.nctl[1]
            cp_xyz[:,i,j,k] = rotMat*cp_xyz[:,i,j,k]
          end
        end
      end


    end
    
    # Rigid body translation
    xfac = 0.2
    yfac = 0.3
    if mesh.dim == 3
      zfac = 0.4
    else
      zfac = 0.0
    end

    if mesh.dim == 2
      cp_xyz[1, :, :] += xfac
      cp_xyz[2, :, :] += yfac
    else
      cp_xyz[1,:,:,:] += xfac
      cp_xyz[2,:,:,:] += yfac
      cp_xyz[3,:,:,:] += zfac
    end

    # set new control point values
    setControlPoints(map, cp_xyz)

    vertices = zeros(Complex128, mesh.dim, map.numFacePts)
    evalSurface(map, vertices)

    for i=1:map.numFacePts
      if mesh.dim == 2
        verts_exact = rotMat*vertices_orig[:, i] + [xfac; yfac]
      else
        verts_exact = rotMat*vertices_orig[:, i] + [xfac; yfac; zfac]
      end

      @fact norm(verts_exact - vertices[:, i]) --> roughly(0.0, atol=1e-13)
    end

#=
    rotMat2 = rotMat[1:2, 1:2]
    for bc=1:length(vertices)
      verts_bc = vertices[bc]
      verts_orig_bc = vertices_orig[bc]

      for i=1:size(verts_bc, 3)
        for j=1:size(verts_bc, 2)
          if mesh.dim == 2
            verts_exact = rotMat[1:2, 1:2]*verts_orig_bc[:, j, i] + [xfac; yfac]
          else
            verts_exact = rotMat*verts_orig_bc[:, j, i] + [xfac; yfac; zfac]
          end

          @fact norm(verts_exact - verts_bc[:, j, i]) --> roughly(0.0, atol=1e-13)
        end
      end
    end
=#

    # test getControlPoints
    cp_xyz2 = getControlPoints(map)
    fill!(cp_xyz2, 0.0)
    getControlPoints(map, cp_xyz2)

    if mesh.dim == 2
      @fact maximum(abs(cp_xyz2 - map._cp_xyz[1:2, :, :, 1])) --> roughly(0.0, atol=1e-13)
    else
      @fact maximum(abs(cp_xyz2 - map._cp_xyz)) --> roughly(0.0, atol=1e-13)
    end

  end  # end facts block
  return nothing
end


function test_jac2{Tffd}(map::PumiMapping{Tffd}, mesh)

  facts("----- Testing Jacobian-vector product -----") do

    Xs_dot = zeros(Tffd, mesh.dim, map.numFacePts)
    Xs_dot2 = zeros(Xs_dot)

    vertices_orig = zeros(Xs_dot)
    evalSurface(map, vertices_orig)

    if Tffd <: Complex
      # test that map.cs produces the same coordinates as map
      map_cs = map.map_cs
      copy!(map_cs._cp_xyz, map._cp_xyz)
      vertices_cs = zeros(Xs_dot)
      evalSurface(map, vertices_cs)

      @fact maximum(abs(vertices_cs - vertices_orig)) --> roughly(0.0, atol=1e-13)

      for i=1:length(map.cp_xyz)
        map._cp_xyz[i] += 1
      end

      copy!(map_cs._cp_xyz, map._cp_xyz)
      evalSurface(map, vertices_orig)
      evalSurface(map_cs, vertices_cs)
      @fact maximum(abs(vertices_cs - vertices_orig)) --> roughly(0.0, atol=1e-13)
    end



    for i=1:10
      fill!(Xs_dot, 0.0)
      fill!(Xs_dot2, 0.0)

      Xcp = getControlPoints(map)
      _Xcp_dot = getControlPoints(map)
      rand!(_Xcp_dot)
      Xcp_dot = real(_Xcp_dot)


      # finite difference
      h = 1e-7
      for j=1:length(Xcp_dot)
        Xcp[j] += Xcp_dot[j]*h
#        map.cp_xyz[j] += Xcp_dot[j]*h
      end

      setControlPoints(map, Xcp)

      vertices_new = zeros(vertices_orig)
      evalSurface(map, vertices_new)

      for j=1:length(Xcp_dot)
        Xcp[j] -= Xcp_dot[j]*h
#        map.cp_xyz[j] -= Xcp_dot[j]*h
      end

      setControlPoints(map, Xcp)

      for j=1:length(Xs_dot)
        Xs_dot[j] = (vertices_new[j] - vertices_orig[j])/h
      end

      # complex step
      evaldXdControlPointProduct(map, Xcp_dot, Xs_dot2)

      for j=1:length(Xs_dot)
        @fact Xs_dot[j] --> roughly(Xs_dot2[j], atol=1e-7)
      end
    end

  end  # end facts block

  return nothing
end


function test_jac(map, mesh)

  facts("----- Testing Transposed Jacobian-vector product -----") do

    # construct entire jacobian
    cp_xyz = getControlPoints(map)
    jac = zeros(mesh.dim*map.numFacePts, length(cp_xyz))
    jac2 = zeros(jac)

    h = 1e-20
    pert = Complex128(0, h)
    Xs = zeros(Complex128, mesh.dim, map.numFacePts)
    for i=1:length(cp_xyz)
      cp_xyz[i] += pert
      setControlPoints(map, cp_xyz)
      evalSurface(map, Xs)

      for j=1:length(Xs)
        jac[j, i] = imag(Xs[j])/h
      end

      cp_xyz[i] -= pert
      setControlPoints(map, cp_xyz)
    end

      # compute reverse mode jacobian
      Xs_dot = zeros(Complex128, mesh.dim, map.numFacePts)
      B_dot = zeros(cp_xyz)
      for i=1:length(Xs_dot)
        Xs_dot[i] = 1
        fill!(B_dot, 0.0)
        evaldXdControlPointTransposeProduct(map, Xs_dot, B_dot)

        for j=1:length(B_dot)
          jac2[i, j] = real(B_dot[j])
        end
        Xs_dot[i] = 0
      end

      @fact norm(jac - jac2) --> roughly(0.0, atol=1e-13)


    println("\nTesting product")
    # test the product using the quadratic form dXs^T dXs/dXcp dXcp
    dXs = real(rand(size(Xs)))
    dXcp = real(rand(size(cp_xyz)))

    # forward product
    for i=1:length(cp_xyz)
      cp_xyz[i] += pert*dXcp[i]
    end
    setControlPoints(map, cp_xyz)
    fill!(Xs, 0.0)
    evalSurface(map, Xs)
    val1 = sum(dXs .* imag(Xs)/h)  # contract with dXs

    for i=1:length(cp_xyz)
      cp_xyz[i] -= pert*dXcp[i]
    end
    setControlPoints(map, cp_xyz)

    # transposed product
    fill!(B_dot, 0.0)
    evaldXdControlPointTransposeProduct(map, dXs, B_dot)
    val2 = sum(dXcp .* B_dot)

    @fact val1 --> roughly(val2, atol=1e-13)


  


  end  # end facts block

  return nothing
end


function runtests()

  mesh, sbp, opts = getTestMesh(2, Complex128)

  # Free Form deformation parameters
  ndim = 2
  order = [3,3,2]  # Order of B-splines in the 3 directions
  nControlPts = [6,6,2]
  offset = [0.25, 0.25, 0.25]
  full_geom = false
  bc_nums = [1] # = opts["BC1"]
#  @assert length(geom_faces) == 1

  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset,
                           bc_nums)

  test_surface(map, mesh)

  # make new objects in case test_surface modified the old ones
  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_jac2(map, mesh)

  # test with real values map

  mesh, sbp, opts = getTestMesh(2, Float64)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_jac2(map, mesh)


  # make new objects in case test_surface modified the old ones
  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)


  test_jac(map, mesh)

  # test with one boundary condition with several geometric entities
  bc_nums = [2]

  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_surface(map, mesh)

  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
  test_jac(map, mesh)

  # test multiple BCs
  bc_nums = [1, 2]
  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_surface(map, mesh)

  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
  test_jac2(map, mesh)


  mesh, sbp, opts = getTestMesh(2, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
  test_jac(map, mesh)



  println("\nTesting 3D")
  # now test 3D
  mesh, sbp, opts = getTestMesh(3, Complex128)

  # Free Form deformation parameters
  ndim = 3
  order = [3,3,3]  # Order of B-splines in the 3 directions
  nControlPts = [6,6,6]
  offset = [0.25, 0.25, 0.25]
  full_geom = false
  bc_nums = [1] # opts["BC1"]
#  @assert length(geom_faces) == 1

  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_surface(map, mesh)

  mesh, sbp, opts = getTestMesh(3, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
                           
  test_jac(map, mesh)


  # test multiple BCs
  bc_nums = [1, 2]

  mesh, sbp, opts = getTestMesh(3, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)

  test_surface(map, mesh)

  mesh, sbp, opts = getTestMesh(3, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
                           
  test_jac2(map, mesh)


  mesh, sbp, opts = getTestMesh(3, Complex128)
  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)
                           
  test_jac(map, mesh)


  return nothing
end

runtests()

#test_surface()
