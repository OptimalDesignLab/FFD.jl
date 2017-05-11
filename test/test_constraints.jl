# test constraints
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
opts["jac_method"] = 2
opts["run_type"] = 5

# Create PumiMesh and SBP objects
sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = PDESolver.createMeshAndOperator(opts, 1)
geom_faces = opts["BC2"]

# Free Form deformation parameters
ndim = 3
order = [2,4,2]  # Order of B-splines in the 3 directions
nControlPts = [2,4,2]

# code here for creating rectilinear map
# Create Mapping object
map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false,
                            geom_faces=geom_faces)

# Create knot vector
calcKnot(map)

# Create Bounding box
offset = [0., 0., 0.5] # No offset in the X & Y direction
ffd_box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

# Control points
controlPoint(map, ffd_box)

# Populate map.xi
calcParametricMappingNonlinear(map, ffd_box, mesh, geom_faces)

facts("--- Checking Linear-Plane Constraints ---") do

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearPlaneConstraints(map, di)
  @fact numCnstr --> 1*3*4

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearPlaneConstraints!(map, di, cntEq, fncidx)
  @fact fncidx --> numCnstr+1
  for k = 1:numCnstr
    @fact cntEq[k] --> 4
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(sum(cntEq))
  Flow = rand(numCnstr)
  Fupp = rand(numCnstr)
  fncidx, ptr = setLinearPlaneConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                           Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @fact iLfun[ptr] --> k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Corner Constraints ---") do

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearCornerConstraints(map, di)
  @fact numCnstr --> 8

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearCornerConstraints!(map, di, cntEq, fncidx)
  @fact fncidx --> numCnstr+1
  for k = 1:numCnstr
    @fact cntEq[k] --> 4
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(sum(cntEq))
  Flow = rand(numCnstr)
  Fupp = rand(numCnstr)
  fncidx, ptr = setLinearCornerConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                            Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @fact iLfun[ptr] --> k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Stretch Constraints ---") do

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearStretchConstraints(map)
  @fact numCnstr --> (2*4*2 - 1)

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearStretchConstraints!(map, cntEq, fncidx)
  @fact fncidx --> numCnstr+1
  for k = 1:numCnstr
    @fact cntEq[k] --> 2
  end

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(sum(cntEq))
  Flow = rand(numCnstr)
  Fupp = rand(numCnstr)
  fncidx, ptr = setLinearStretchConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                             Fupp, fncidx, ptr)

  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @fact iLfun[ptr] --> k
      ptr += 1
    end
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Root Constraints ---") do

  di = 2

  # test function for counting number of constraints
  numCnstr = numLinearRootConstraints(map)
  @fact numCnstr --> 5

  # test function for counting variables in each constraint
  cntEq = zeros(Int, numCnstr)
  fncidx = 1
  fncidx = countVarsLinearRootConstraints!(map, cntEq, fncidx)
  @fact fncidx --> numCnstr+1
  for k = 1:3
    @fact cntEq[k] --> 1
  end
  @fact cntEq[4] --> 2
  @fact cntEq[5] --> 2

  # test function for filling sparse constraint Jacobian
  fncidx = 1
  ptr = 1
  iLfun = zeros(Int, sum(cntEq))
  jLvar = zeros(Int, sum(cntEq))
  LinG = zeros(sum(cntEq))
  Flow = rand(numCnstr)
  Fupp = rand(numCnstr)
  fncidx, ptr = setLinearRootConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                          Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  ptr = 1
  for k = 1:numCnstr
    for i = 1:cntEq[k]
      @fact iLfun[ptr] --> k
      ptr += 1
    end
  end
  for k = 1:3
    @fact LinG[k] --> 1.0
  end
  @fact LinG[4] --> 1.0
  @fact LinG[5] --> -1.0
  @fact LinG[6] --> 1.0
  @fact LinG[7] --> -1.0
end

MPI.Finalize()
