# Type Mapping
@doc """
### Mapping

A type for creating mapping objects needed for creating

"""->

type Mapping
  nctl::Array(Int, 3)   # Number of control points in each of the 3 dimensions
  jkmmax::Array(Int, 3) # Number of nodes in each direction
  order::Array(Int, 3)  # Order of B-spline in each direction

  xi = AbstractArray{AbstractFloat, 4}     # Coordinate values
  cp_xyz = AbstractArray{AbstractFloat, 4} # Cartesian coordinates of control points
  edge_knot = AbstractArray{AbstractFloat, 3}  # edge knot vectors
  edge_param = AbstractArray{AbstractFloat, 3} # edge parameters

  # Working arrays
  aj = AbstractArray{AbstractFloat, 3}
  dl = AbstractArray{AbstractFloat, 2}
  dr = AbstractArray{AbstractFloat, 2}
  knot = AbstractArray{AbstractFloat, 2}
  work = AbstractArray{AbstractFloat, 4}

  function Mapping(ncpts, nnodes, k)

    # Define max_wrk = number of work elements at each control point
    const max_work = 2*6  # 2*n_variables

    nctl[:] = ncpts[:]    # Set number of control points
    jkmmax[:] = nnodes[:] # Set map refinement level in each coordinate direction
    order[:] = k[:]       # Set the order of B-splines

    # Assertion statements to prevent errors
    for i = 1:3
      @assert nctl[i] > 0
      @assert jkmmax[i] > 0
      @assert order[i] > 0
    end

    # Allocate and initialize mapping arrays
    max_order = maximum(order)  # Highest order among 3 dimensions
    max_knot = max_order + maximum(nctl) # Maximum number of knots among 3 dimensions

    xi = zeros(jkmmax[1], jkmmax[2], jkmmax[3], 3)
    cp_xyz = zeros(nctl[1], nctl[2], nctl[3], 3)
    edge_knot = zeros(max_knot, 4, 3)
    knot = zeros(max_knot, 3)
    edge_param = zeros(2,4,3)
    aj = zeros(3, max_order, 3)
    dl = zeros(max_order-1, 3)
    dr = zeros(max_order-1, 3)
    work = zeros(ntcl[1], nctl[2], nctl[3], max_work)

    new(nctl, jkmmax, order, xi, cp_xyz, edge_knot, edge_param, aj, dl, dr,
        knot, work)

  end  # End constructor

end  # End Mapping
