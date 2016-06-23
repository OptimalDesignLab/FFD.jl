# test components


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
map = LinearMapping(ndim, order, nControlPts, nnodes)

# Create BoundingBox object
offset = [0.5,0.5,0.5]  # offset for the bounding box
geom_bounds = [1. 1. 1.;3. 3. 3.]
box = BoundingBox(ndim, geom_bounds, offset)

facts("--- Checking Mapping object ---") do

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

facts("--- Checking BoundingBox ---") do

  @fact box.ndim --> 3
  @fact box.origin --> [0.5, 0.5, 0.5]
  @fact box.unitVector --> [1. 0. 0.;0. 1. 0.;0. 0. 1.]
  @fact box.geom_coord --> [1. 1. 1.;3. 3. 3.]
  @fact box.offset --> [0.5,0.5,0.5]
  @fact box.box_bound --> [0.5 0.5 0.5;3.5 3.5 3.5]

end  # End facts("--- Checking BoundingBox ---")

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

facts("--- Checking FFD Volume Evaluation with Linear Mapping---") do

  Vol = zeros(nodes_xyz)
  evalVolume(map, Vol)
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

end  # facts("--- Checking FFD Volume Evaluation ---")

facts("--- Checking Nonlinear Mapping ---") do

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

  context("Check nonlinearMap") do

    X = ones(AbstractFloat, 3)
    pX = zeros(X)

    nonlinearMap(map, box, X, pX)
    for i = 1:3
      @fact pX[i] --> roughly(0.16666666666666666, atol = 1e-15)
    end

  end  # End context ("Check nonlinearMap")

end # End facts("--- Checking Linear Mapping ---")
