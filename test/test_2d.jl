
push!(LOAD_PATH, "../src/")

using FreeFormDeformation
using PumiInterface
using PdePumiInterface
using SummationByParts
using ODLCommonTools
using FactCheck


function getTestMesh{Tmsh}(shape_type::Integer, ::Type{Tmsh})

  degree = 1

  if shape_type == 2
    sbp = getTriSBPOmega(degree=degree)
  elseif shape_type == 3
    sbp = getTriSBPGamma(degree=degree)
  end

  sbpface = TriFace{Float64}(degree, sbp.cub, sbp.vtx)

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

  mesh = PumiMeshDG2(Tmsh, sbp, opts, sbpface, dofpernode=1,
                    shape_type=shape_type)

  return mesh, sbp, opts
end


function test_surface()


  facts("----- Testing FFD Surface Evaluation -----") do
    mesh, sbp, opts = getTestMesh(2, Complex128)


    # Free Form deformation parameters
    ndim = 2
    order = [2,2,2]  # Order of B-splines in the 3 directions
    nControlPts = [2,2,2]
    offset = [0.25, 0.25, 0.25]
    full_geom = false
    geom_faces = opts["BC1"]
    @assert length(geom_faces) == 1

    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, full_geom,
                             geom_faces)

    vertices_orig = evalSurface(map, mesh)

    # Rigid body rotation
    theta = -20*pi/180  # Rotate wall coordinates by 10 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
  
    cp_orig = copy(map.cp_xyz)
    # Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # Rigid body translation
    xfac = 0.2
    yfac = 0.3
    map.cp_xyz[1,:,:,:] += xfac
    map.cp_xyz[2,:,:,:] += yfac

    vertices = evalSurface(map, mesh)

    rotMat2 = rotMat[1:2, 1:2]
    for bc=1:length(vertices)
      verts_bc = vertices[bc]
      verts_orig_bc = vertices_orig[bc]

      for i=1:size(verts_bc, 3)
        for j=1:size(verts_bc, 2)
          verts_exact = rotMat2*verts_orig_bc[:, j, i] + [xfac; yfac]
          @fact norm(verts_exact - verts_bc[:, j, i]) --> roughly(0.0, atol=1e-13)
        end
      end
    end
  end  # end facts block

  return nothing
end


function test_jac()

  facts("----- Testing Transposed Jacobian-vector product -----") do
    mesh, sbp, opts = getTestMesh(2, Complex128)


    # Free Form deformation parameters
    ndim = 2
    order = [2,2,2]  # Order of B-splines in the 3 directions
    nControlPts = [2,2,2]
    offset = [0.25, 0.25, 0.25]
    full_geom = false
    geom_faces = opts["BC1"]
    @assert length(geom_faces) == 1

    map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, full_geom,
                             geom_faces)


    # compute entire jacbian
    Xs_orig = evalSurface(map, mesh)
    nXs = length(Xs_orig[1])
    nCP = length(map.cp_xyz)
    jac = zeros(nXs, nCP)
    jac2 = zeros(jac)

    # compute complex step jacobian
    h = 1e-20
    pert = Complex128(0, h)
    for i=1:nCP
      map.cp_xyz[i] += pert
      Xs = evalSurface(map, mesh)

      for j=1:nXs
        jac[j, i] = imag(Xs[1][j])/h
      end

      map.cp_xyz[i] -= pert
    end

    # compute reverse mode jacobian
    Xs_dot = Array(Array{Complex128, 3}, 1)
    Xs_dot[1] = zeros(size(Xs_orig[1]))
    B_dot = zeros(map.cp_xyz)
    for i=1:nXs
      Xs_dot[1][i] = 1
      fill!(B_dot, 0.0)
      FreeFormDeformation.evaldXdControlPointTransposeProduct(map, mesh, Xs_dot, B_dot)

      for j=1:nCP
        jac2[i, j] = real(B_dot[j])
      end
      Xs_dot[1][i] = 0
    end

    @fact norm(jac - jac2) --> roughly(0.0, atol=1e-13)
  end

  return nothing
end


test_surface()
test_jac()
