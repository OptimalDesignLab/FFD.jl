# Type Mapping

type Mapping
  nctl::Array(Int, 3)  # Number of control points in each of the 3 dimensions
  jkmax::Array(Int, 3) # Number of nodes in each direction
  order::Array(Int, 3) # Order of B-spline in each direction

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

end  # End Mapping
