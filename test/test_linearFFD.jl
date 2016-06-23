# test components

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
    calcdXdxi(map, xi, jderiv, dX)
    @fact dX[1] --> roughly(3.0, atol = 1e-15)
    @fact dX[2] --> roughly(0.0, atol = 1e-15)
    @fact dX[3] --> roughly(0.0, atol = 1e-15)

  end  # End context("Checking Volume Derivative")

  context("Checking nonlinearMap") do

    X = ones(AbstractFloat, 3)
    pX = zeros(X)

    nonlinearMap(map, box, X, pX)
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
    evalVolumePoint(map, xi, xyz)
    for i = 1:map.ndim
      @fact xyz[i] --> roughly(2.0, atol = 1e-15)
    end

  end  # End context("Checking single poitn evaluation")

  context("Checking multiple point evaluation") do
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

  end # End context("Checking multiple point evaluation")

end  # facts("--- Checking FFD Volume Evaluation ---")
