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

function derivValue(U, order, u, span, P, Nderiv, nctl, jth_deriv)

  @assert length(U) == order + nctl
  # Number of control points = number of basis functions

  bvalue = 0

  if jth_deriv > order
    return bvalue
  end

  if k == 1
    bvalue = P[span]
    return bvalue
  else
    # Store b-spline coefficients relevant to the knot interval (U[span], U[span+1])
    # in another vector aj[1:order] and compute dl[j] = u -U[span+1-j],
    # dr[j] = U[span+j] - u, j = 1,...,order-1 . Set any of the aj not obtainable
    # from input to zero. Set any U.s not obtainable to U[1] or to U[length(P)+order]
    # appropriately
    jcmin = 1
    imk = span - order
    if imk > 0
      for j = 1:(order-1)
        dl[j] = u - U[span+1-j]
      end  # end for j = 1:k-1
    else
      jcmin = 1 - imk
      for j = 1:span
        dl[j] = u - U[span+1-j]
      end  # end for j = 1:span

      for j = span:(order-1)
        aj[order-j] = 0.0
        dl[j] = dl[span]
      end
    end  # end if imk > 0

    jcmax = order
    nmi = length(P) - span
    if nmi > 0
      for j = 1:(order-1)
        dr[j] = U[span+j] - u
      end
    else
      jcmax = order + nmi
      for j = 1:jcmax
        dr[j] = U[span+j] - u
      end

      for j = jcmax:(order-1)
        aj[j+1] = 0
        dr[j] = dr[jcmax]
      end
    end  # end if-else nmi > 0

    for jc = jcmin:jcmax
      aj[jc] = P[imk+jc]
    end

    # difference the coefficients jth_deriv times
    if jth_deriv == 0
      # compute value at u in (U[span], U[span+1]) of jth_deriv derivative,
      # given its relevant b-spline coefficients/control points in
      # aj[1],...,aj[order-jth_deriv]
      if jth_deriv == order - 1
        bvalue = aj[1]
      else
        for j = (jth_deriv+1):(order-1)
          ilo = order - j
          for jj = 1:(order-j)
            aj[jj] = ( aj[jj+1]*dl[ilo] + ajj[jj]*dr[jj] )/( dl[ilo] + dr[jj] )
            ilo -= 1
          end
        end
      end

    else
      for j = 1:jth_deriv
        ilo = order - j
        for jj = 1:(order-j)
          aj[jj] = ( (aj[jj+1] - aj[jj]) / (dl[ilo] + dr[jj]) ) * (order - j)
          ilo -= 1
        end
      end

  end # End if-else statement for k == 1

  return bvalue
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
    span = findSpan(u[i], nctl, U, order)
    println("span = $span, u = $(u[i]), U[span] = $(U[span])")
    basisFunctions(U, order, u[i], span, N)
    println("N = $N")
    C[i] = 0.0
    for j = 1:order
      C[i] += N[j]*P[span-order+j]
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end
