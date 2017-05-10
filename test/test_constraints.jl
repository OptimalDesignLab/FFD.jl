# test constraints

facts("--- Checking Linear-Plane Constraints ---") do

  # Free Form deformation parameters
  ndim = 3
  order = [2,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,4,2]

  # code here for creating rectilinear map
  
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
  Fupp = rand(numCntrs)
  fncidx, ptr = setLinearPlaneConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                           Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  for k = 1:numCnstr
    @fact iLfnc[k] --> k
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Corner Constraints ---") do

  # Free Form deformation parameters
  ndim = 3
  order = [2,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,4,2]

  # code here for creating rectilinear map
  
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
  Fupp = rand(numCntrs)
  fncidx, ptr = setLinearCornerConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                            Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  for k = 1:numCnstr
    @fact iLfnc[k] --> k
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Stretch Constraints ---") do

  # Free Form deformation parameters
  ndim = 3
  order = [2,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,4,2]

  # code here for creating rectilinear map
  
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
  Fupp = rand(numCntrs)
  fncidx, ptr = setLinearStretchConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                             Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  @fact Flow --> roughly(zeros(numCnstr), 1e-15)
  @fact Fupp --> roughly(zeros(numCnstr), 1e-15)
  for k = 1:numCnstr
    @fact iLfnc[k] --> k
  end
  # !!!!!!!!!! should also include checks on jLvar and LinG, eventually
end

facts("--- Checking Linear-Root Constraints ---") do

  # Free Form deformation parameters
  ndim = 3
  order = [2,4,2]  # Order of B-splines in the 3 directions
  nControlPts = [2,4,2]

  # code here for creating rectilinear map
  
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
  Fupp = rand(numCntrs)
  fncidx, ptr = setLinearRootConstraints!(map, di, iLfun, jLvar, LinG, Flow,
                                          Fupp, fncidx, ptr)
  @fact fncidx --> numCnstr+1
  @fact ptr --> sum(cntEq)+1
  for k = 1:numCnstr
    @fact iLfnc[k] --> k
  end
  for k = 1:3
    @fact LinG[k] --> 1.0
  end
  @fact LinG[4] --> 1.0
  @fact LinG[5] --> -1.0
  @fact LinG[6] --> 1.0
  @fact LinG[7] --> -1.0
end
