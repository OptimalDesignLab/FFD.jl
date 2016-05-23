# Functions based on de boor book

function basisFunctions(T, jhigh, index, x, left, N)

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
