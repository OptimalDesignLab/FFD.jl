# B-Spline Basis function for Bezier Bernsetirn splines
@doc """
### basisFunctions

Computes the basis functions needed to compute B-Splines

**Inputs**

*  `i` : Knot span index
*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `p` : Degree of B-spline basis function (Order p+1)
*  `U` : Knot vector
*  `N` : Array of basis functions

**Outputs**

*  None

SOURCE: The NURBS book 2nd Edition, Algorithm A2.2

"""->

function basisFunctions(i, u, p, U, N)

  # Initialize variables
  left = Array(Float64, p+1)
  right = Array(Float64, p+1)
  saved = 0.0

  N[1] = 1.0

  for j = 2:p+1
    left[j] = u - U[i+1-j]
    right[j] = U[i+j] - u
    saved = 0.0

    for r = 1:j-1
      temp = N[r]/(right[r+1] + left[j-r])
      N[r+1] = saved + right[r+1]*temp
      saved = left[r]*temp
    end

    N[j] = saved
  end
  println("\nN = $N\n")
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
*  `p` : Degree of B-spline basis function, (Order p+1)
*  `P` : Array of control points
*  `C` : Resulting curve values

**Outputs**

*  None

SOURCE: The NURBS book 2nd Edition, Algorithm A3.1
      : Gaetan's pyspline/src/eval_curve
"""->

function evalCurve(u, U, p, P, C)

  nctl = length(P)
  N = Array(Float64, p+1) # Array of basis functions 1D
  for i = 1:length(u)
    span = findSpan(u[i], nctl, U, p)
    println("\nspan = $span\n")
    basisFunctions(span, u[i], p, U, N)
    C[i] = 0.0
    for j = 1:p+1
      C[i] += N[i]*P[span-p+i]
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end
