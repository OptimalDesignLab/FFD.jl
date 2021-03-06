# test constraints
# MPI Declarations
comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

# opts = PdePumiInterface.get_defaults()
# # 2D mesh
# opts["order"] = 1
# opts["dimensions"] = 2
# opts["use_DG"] = true
# opts["operator_type"] = "SBPOmega"
# opts["dmg_name"] = "../src/mesh_files/2D_Airfoil.dmg"
# opts["smb_name"] = "../src/mesh_files/2D_Airfoil.smb"
# opts["numBC"] = 2
#
# # For 2DAirfoil
# opts["BC1"] = [8,11,14,17]
# opts["BC1_name"] = "FarField"
# opts["BC2"] = [5]
# opts["BC2_name"] = "Airfoil"
#
# opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
# opts["jac_type"] = 2
# opts["jac_method"] = 2
# opts["run_type"] = 5
#
# # Create PumiMesh and SBP objects
# sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = PDESolver.createMeshAndOperator(opts, 1)

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

#mesh = PumiMeshDG2{Tmsh}(opts["dmg_name"], opts["smb_name"], opts["order"],
#                         sbp, opts, sbpface; dofpernode=dofpernode,
#                         coloring_distance=opts["coloring_distance"],
#                         shape_type=shape_type)

#geom_faces = opts["BC2"]
bc_nums = [2]
# Free Form deformation parameters
ndim = 3
order = [2,4,2]  # Order of B-splines in the 3 directions
nControlPts = [2,4,2]

# code here for creating rectilinear map
# Create Mapping object
map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh,
                            bc_nums=bc_nums)

# Create knot vector
calcKnot(map)

# Create Bounding box
offset = [0., 0., 0.5] # No offset in the X & Y direction
ffd_box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

# Control points
controlPoint(map, ffd_box)

# Populate map.xi
calcParametricMappingNonlinear(map, ffd_box, mesh, bc_nums)

@testset "--- Checking Linear-Plane Constraints ---" begin

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearPlaneConstraints(map, di)
  @test ( numCnstr )== 1*3*4

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearPlaneConstraints!(map, di, cntEq, fncidx)
  @test ( fncidx )== numCnstr+1
  for k = 1:numCnstr
    @test ( cntEq[k] )== 4
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(Tmsh, sum(cntEq))
  Flow = rand(numCnstr) + 0im
  Fupp = rand(numCnstr) + 0im
  fncidx, ptr = setLinearPlaneConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                           Fupp, fncidx, ptr)
  @test ( fncidx )== numCnstr+1
  @test ( ptr )== sum(cntEq)+1
  @test isapprox( Flow, zeros(numCnstr), 1e-15) 
  @test isapprox( Fupp, zeros(numCnstr), 1e-15) 
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @test ( iLfun[ptr] )== k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

@testset "--- Checking Linear-Corner Constraints ---" begin

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearCornerConstraints(map, di)
  @test ( numCnstr )== 8

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearCornerConstraints!(map, di, cntEq, fncidx)
  @test ( fncidx )== numCnstr+1
  for k = 1:numCnstr
    @test ( cntEq[k] )== 4
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(Tmsh, sum(cntEq))
  Flow = rand(numCnstr) + 0im
  Fupp = rand(numCnstr) + 0im
  fncidx, ptr = setLinearCornerConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                            Fupp, fncidx, ptr)
  @test ( fncidx )== numCnstr+1
  @test ( ptr )== sum(cntEq)+1
  @test isapprox( Flow, zeros(numCnstr), 1e-15) 
  @test isapprox( Fupp, zeros(numCnstr), 1e-15) 
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @test ( iLfun[ptr] )== k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

@testset "--- Checking Linear-Stretch Constraints ---" begin

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearStretchConstraints(map)
  @test ( numCnstr )== (2*4*2 - 1)

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearStretchConstraints!(map, cntEq, fncidx)
  @test ( fncidx )== numCnstr+1
  for k = 1:numCnstr
    @test ( cntEq[k] )== 2
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(Tmsh, sum(cntEq))
  Flow = rand(numCnstr) + 0im
  Fupp = rand(numCnstr) + 0im
  fncidx, ptr = setLinearStretchConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                             Fupp, fncidx, ptr)

  @test ( fncidx )== numCnstr+1
  @test ( ptr )== sum(cntEq)+1
  @test isapprox( Flow, zeros(numCnstr), 1e-15) 
  @test isapprox( Fupp, zeros(numCnstr), 1e-15) 
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @test ( iLfun[ptr] )== k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

@testset "--- Checking Linear-Root Constraints ---" begin

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearRootConstraints(map)
  @test ( numCnstr )== 5

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearRootConstraints!(map, cntEq, fncidx)
  @test ( fncidx )== numCnstr+1
  for k = 1:3
    @test ( cntEq[k] )== 1
  end
  @test ( cntEq[4] )== 2
  @test ( cntEq[5] )== 2

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(Tmsh, sum(cntEq))
  Flow = rand(numCnstr) + 0im
  Fupp = rand(numCnstr) + 0im
  fncidx, ptr = setLinearRootConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                          Fupp, fncidx, ptr)
  @test ( fncidx )== numCnstr+1
  @test ( ptr )== sum(cntEq)+1
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @test ( iLfun[ptr] )== k
      ptr += 1
    end
  end
  for k = 1:3
    @test ( LinG[k] )== 1.0
  end
  @test ( LinG[4] )== 1.0
  @test ( LinG[5] )== -1.0
  @test ( LinG[6] )== 1.0
  @test ( LinG[7] )== -1.0
end
