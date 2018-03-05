
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


function test_jac()

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

  writedlm("jac.dat", jac)
  writedlm("jac2.dat", jac2)
  println("jac = \n", jac)
  println("jac2 = \n", jac2)
  println("diff = \n", jac - jac2)
  println("diffnorm = \n", norm(jac - jac2))

  return nothing
end


test_jac()
