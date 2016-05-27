# B-Spline Basis function for Bezier Bernsetirn splines
@doc """
### basisFunctions

Computes the basis functions needed to compute B-Splines. The calculation starts
from scratch and the entire triangular array of B-spline values of orders
1,2,...,jhigh is generated order by order (column by column).

**Inputs**

*  `U`     : Knot vector
*  `order` : Order of B-spline
*  `u`     : Location where basis functions are to be evaluated
*  `span`  : Knot span index
*  `N`     : Bassi Functions

**Outputs**

*  None

SOURCE: "A practical guide to splines" bt C. de Boor Page 112 BSPLVB

"""->
function basisFunctions(U, order, u, span, N)

  jmax = 20 # arbitrarily imposed
  delta_l = zeros(Float64, jmax)
  delta_r = zeros(Float64, jmax)

  N[1] = 1.0

  for j = 1:(order-1)
    delta_r[j] = U[span+j] - u
    delta_l[j] = u - U[span+1-j]
    saved = 0.0
    for i = 1:j
      term = N[i]/(delta_r[i] + delta_l[j+1-i])
      N[i] = saved + delta_r[i]*term
      saved = delta_l[j+1-i]*term
    end  # end for i = 1:r
    N[j+1] = saved
  end  # end for j = 1:(order-1)

  return nothing
end

@doc """
### derivBasisFunctions
Computes the derivatives of the basis functions at a particular point u.

**Inputs**

*  `map` : Object of mapping type
*  `N`   : Basis functions
*  `Nderiv` : Derivatives of the basis functions at u
*  `U`   : Knot vector
*  `di`  : Direction in which the derivative is to be evaluated
*  `u`   : Location where the derivatives are to be calculated
*  `span`: Knot span index of u

**Outputs**

*  None

"""->

function derivBasisFunctions(map, N, Nderiv, di, u, span)

  order = map.order[di]
  N[1] = 1.0
  Nderiv[1] = 0.0

  dr = view(map.dr, :, :)
  dl = view(map.dr, :, :)

  if order > 1
    for k = 1:order-1
      kp1 = k + 1
      dr[k,di] = map.knot[span+k,di] - u
      dl[k,di] = u - map.knot[span+1-k, di]
      saved = 0.0
      dsaved = 0.0
      for i = 1:k
        temp = 1/( dr[i, di] + dl[k+1-i, di] )
        dtemp = Nderiv[i]*temp
        temp = N[i]*temp

        N[i] = saved + dr[i,di]*temp
        Nderiv[i] = dsaved - temp + dr[i,di]*dtemp

        saved = dl[k+1-i,di]*temp
        dsaved = temp + dl[k+1-i, di]*dtemp
      end  # end for i = 1:k
      N[k+1] = saved
      Nderiv[k+1] = dsaved
    end  # end for k = 1:order-1
  end  # end if order > 1

  return nothing
end

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

  if u >= U[n+1]
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
