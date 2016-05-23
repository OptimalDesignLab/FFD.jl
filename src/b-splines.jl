# B-Spline Basis function for Bezier Bernsetirn splines
@doc """
### basisFunctions

Computes the basis functions needed to compute B-Splines

**Inputs**

*  `i` : Knot span index
*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `p` : Degree of B-spline basis function (Order p+1)
*  `U` : Knot vector
*  `N` : Array of basis functions or the B-spline coefficients

**Outputs**

*  None

SOURCE: "A practical guide to splines" bt C. de Boor Page 112 BSPLVB

"""->

function basisFunctions(i, u, p, U, N)

  # Initialize variables
  jmax = 20
  delta_l = Array(Float64, jmax)
  delta_r = Array(Float64, jmax)
  saved = 0.0
  # jhigh = p+1
  # println("jhigh = $jhigh")

  N[1] = 1.0
  for j = 1:p
    delta_l[j] = U[i+j] - u
    delta_r[j] = u - U[i+1-j]
    saved = 0.0
    for r = 1:j
      term = N[r]/(delta_r[r] + delta_l[j+1-r])
      N[r] = saved + delta_r[r]*term
      saved = delta_l[j+1-r]*term
    end
    N[j+1] = saved
  end
  # println("N = $N")
  return nothing
end  # End function basisFunctions

@doc """
### findSpan

Determines the knot span index, i

**Inputs**

*  `n` : Number of control points
*  `p` : Degree of B-spline basis function (Order p+1)
*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `U` : Knot vector

**Outputs**

*  `span` : Knot span index

SOURCE: The NURBS book 2nd Edition, Algorithm A2.1

"""->

function findSpan(u, n, U, p)

  if u == U[n+1]
    return n  # Special case
  end

  low = p
  high = n+1
  mid = div(low + high, 2)
  # Do a binary search
  while u < U[mid] || u >= U[mid+1]
    if u < U[mid]
      high = mid
    else
      low = mid
    end  # End if
    mid = div(low+high, 2)
  end  # End while

  return mid
end  # End function findSpan

@doc """
### evalCurve

Determine the value of curve at certain points

**Inputs**

*  `u` : Array of coordinates where the curve needs to be computed
*  `U` : Knot vector
*  `order` : order of B-spline basis functio, (order = p+1)
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
  N = Array(Float64, p+1) # Array of basis functions 1D
  for i = 1:length(u)
    span = findSpan(u[i], nctl, U, p)
    println("span = $span, u = $(u[i])")
    basisFunctions(span, u[i], p, U, N)
    C[i] = 0.0
    for j = 1:p+1
      println(span-p+j-2)
      C[i] += N[j]*P[span-p+j-2] # -2 because compared to algorithm A3.1 this
                                 # this uses 1 based indexeng. hence (span-1) &
                                 # (j-1) when compared to A3.1
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end
