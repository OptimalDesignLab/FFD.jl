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
function basisFunctions(T, jhigh, x, left, N)

  # T = knot vector
  # jhigh = Max number of bases
  # index = integer that determines the order. jout = max(jhigh, (j+1)*(index-1))
  #         of the bsplines whose values at x are to be returned. Index is used
  #         to avoid recalculations when several columns of the triangular Array
  #         of B-spline values are needed
  # x = location where bases are evaluated
  # left = knot span index
  # N = bases
  # The calculation starts from scratch and the entire triangular array of B-spline
  # values of orders 1,2,...,jhigh is generated order by order (column by column)

  jmax = 20 # arbitrarily imposed
  delta_l = zeros(Float64, jmax)
  delta_r = zeros(Float64, jmax)

  jout = length(T) - left

  N[1] = 1.0

  for j = 1:(jhigh-1)
    delta_r[j] = T[left+j] - x
    delta_l[j] = x - T[left+1-j]
    saved = 0.0
    for i = 1:j
      term = N[i]/(delta_r[i] + delta_l[j+1-i])
      N[i] = saved + delta_r[i]*term
      saved = delta_l[j+1-i]*term
    end  # end for i = 1:r
    N[j+1] = saved
  end  # end for j = 1:jhigh
  # println("N = $N")
  return nothing
end

#=
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
    println("delta_l = $delta_l")
    println("delta_r = $delta_r")
    saved = 0.0
    for r = 1:j
      term = N[r]/(delta_r[r] + delta_l[j+1-r])
      println("term = $term")
      N[r] = saved + delta_r[r]*term
      saved = delta_l[j+1-r]*term
      println("saved = $saved")
      println("In basis function, N[r] = $(N[r]), N = $N")
    end
    N[j+1] = saved
    println("\n\n")
  end
  # println("N = $N")
  return nothing
end  # End function basisFunctions
=#
@doc """
### findSpan

Determines the knot span index, i

**Inputs**

*  `n` : Number of control points
*  `k` : Order of B-spline basis function
*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `U` : Knot vector

**Outputs**

*  `span` : Knot span index

SOURCE: The NURBS book 2nd Edition, Algorithm A2.1

"""->

function findSpan(u, n, U, k)

  if u == U[n+1]
    return n  # Special case when u = last term of knot vector
  elseif u < U[k]
    return k  # When u lies at the starting point of the curve
  end

  low = k-1
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
    span = findSpan(u[i], nctl, U, order)
    println("span = $span, u = $(u[i]), U[span] = $(U[span])")
    # basisFunctions(span, u[i], p, U, N)
    basisFunctions(U, order, u[i], span, N)
    println("N = $N")
    C[i] = 0.0
    for j = 1:order
      # println(span - order + j)
      C[i] += N[j]*P[span-order+j] # -2 because compared to algorithm A3.1 this
                                 # this uses 1 based indexeng. hence (span-1) &
                                 # (j-1) when compared to A3.1
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end
