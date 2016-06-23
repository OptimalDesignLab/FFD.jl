# Evaluations

@doc """
### evalCurve

Determine the value of curve at certain points

**Inputs**

*  `u` : Array of coordinates where the curve needs to be computed
*  `U` : Knot vector
*  `order` : order of B-spline basis functions, (order = p+1)
*  `P` : Array of control points
*  `C` : Resulting curve values

**Outputs**

*  None

SOURCE: The NURBS book 2nd Edition, Algorithm A3.1
      : Gaetan's pyspline/src/eval_curve
"""->

function evalCurve(u, U, order, P, C)

  @assert length(P) + order == length(U)
  p = order - 1  # Degree of B-spline basis function
  nctl = length(P)
  N = Array(Float64, order) # Array of basis functions 1D

  for i = 1:length(u)
    span = findSpan(u[i], U, order, nctl)
    basisFunctions(U, order, u[i], span, N)
    C[i] = 0.0
    for j = 1:order
      C[i] += N[j]*P[span-order+j]
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end

@doc """
### evalVolume

Determine the the coordinates of all points in a 3D volume using B-splines. The
symbol convention used in the function is from the book

"The NURBS book 2nd Edition"

**Arguments**

*  `map` : Object of mapping type
*  `Vol` : (x,y,z) coordinates of the embedded volume within the contol points

"""->

function evalVolume(map, Vol)

  fill!(Vol, 0.0) # Zero out all entries of Vol

  for k = 1:map.numnodes[3]
    for j = 1:map.numnodes[2]
      for i = 1:map.numnodes[1]
        xyz = view(Vol, i,j,k,:)
        evalVolumePoint(map, map.xi[i,j,k,:], xyz)
      end  # End for i = 1:map.numnodes[3]
    end    # End for j = 1:map.numnodes[2]
  end      # End for k = 1:map.numnodes[1]

  return nothing
end

@doc """
### evalVolumePoint

Computes a pont in the FFD volume. The symbol convention used in this function
is from "The NURBS book 2nd Edition"

**Arguments**

*  `map` : Object of Mapping type
*  `xi`  : Parametric FFD coordinates (referred to as the (s,t,u) coordinate
           systems within FFD functions and as (u,v,w) in this function)
*  `xyz` : (x,y,z) coordinates of the parametric point in the FFD volume
"""->

function evalVolumePoint(map, xi, xyz)

  Nu = zeros(map.order[1])
  Nv = zeros(map.order[2])
  Nw = zeros(map.order[3])

  # Work with u
  span = findSpan(xi[1], map.edge_knot[1], map.order[1], map.nctl[1])
  basisFunctions(map.edge_knot[1], map.order[1], xi[1], span, Nu)
  startu = span - map.order[1]

  # Work with v
  span = findSpan(xi[2], map.edge_knot[2], map.order[2], map.nctl[2])
  basisFunctions(map.edge_knot[2], map.order[2], xi[2], span, Nv)
  startv = span - map.order[2]

  # Work with w
  span = findSpan(xi[3], map.edge_knot[3], map.order[3], map.nctl[3])
  basisFunctions(map.edge_knot[3], map.order[3], xi[3], span, Nw)
  startw = span - map.order[3]

  for ii = 1:map.order[1]
    for jj = 1:map.order[2]
      for kk = 1:map.order[3]
        for idim = 1:map.ndim
          xyz[idim] += Nu[ii]*Nv[jj]*Nw[kk]*
                       map.cp_xyz[startu+ii, startv+jj, startw+kk, idim]
        end
      end  # End for kk = 1:map.order[3]
    end    # End for jj = 1:map.order[2]
  end      # End for ii = 1:map.order[1]

  return nothing
end
