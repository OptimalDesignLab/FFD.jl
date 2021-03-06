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
*  `N`     : Basis Functions

**Outputs**

*  None

SOURCE: "A practical guide to splines" bt C. de Boor Page 112 BSPLVB

"""->
function basisFunctions(U::AbstractVector{T}, order, u::T2, span, N::AbstractVector{T}) where {T, T2}

  # (U, order, u, span, N)
  #order = map.order[di]
  # U = map.edge_knot[di]

  Tffd = promote_type(T, T2)
  jmax = 20 # arbitrarily imposed
  delta_l = zeros(Tffd, jmax)
  delta_r = zeros(Tffd, jmax)

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

Computes the basis functions and its first & second (method 2) derivative at a
particular point u.

**Inputs**

*  `U`     : Knot vector
*  `order` : Order of B-spline
*  `u`     : Location where basis functions are to be evaluated
*  `span`  : Knot span index
*  `N`     : Basis Functions
*  `Nderiv`: Derivative of the basis function

**Outputs**

*  None

"""->

function derivBasisFunctions(u::T, U::AbstractVector{T2},
                             order, span, N, Nderiv) where {T, T2}

  Tffd = promote_type(T, T2)

  # order = map.order[di]
  N[1] = 1.0
  Nderiv[1] = 0.0

  # dr =sview(map.dr, :, :)
  # dl =sview(map.dl, :, :)

  dr = zeros(Tffd, order-1)
  dl = zeros(Tffd, order-1)

  if order > 1
    for k = 1:order-1
      kp1 = k + 1
      dr[k] = U[span+k] - u
      dl[k] = u - U[span+1-k]
      saved = 0.0
      dsaved = 0.0
      for i = 1:k
        temp = 1/( dr[i] + dl[k+1-i] )
        dtemp = Nderiv[i]*temp
        temp = N[i]*temp

        N[i] = saved + dr[i]*temp
        Nderiv[i] = dsaved - temp + dr[i]*dtemp

        saved = dl[k+1-i]*temp
        dsaved = temp + dl[k+1-i]*dtemp
      end  # end for i = 1:k
      N[k+1] = saved
      Nderiv[k+1] = dsaved
    end  # end for k = 1:order-1
  end  # end if order > 1

  return nothing
end

function derivBasisFunctions(map, N, Nderiv, N2deriv, di, u, span)

  order = map.order[di]
  N[1] = 1.0
  Nderiv[1] = 0.0
  N2deriv[1] = 0.0

  dr =sview(map.dr, :, :)
  dl =sview(map.dl, :, :)

  if order > 1
    for k = 1:order-1
      dr[k,di] = map.knot[span+k,di] - u
      dl[k,di] = u - map.knot[span+1-k, di]
      saved = 0.0
      dsaved = 0.0
      d2saved = 0.0
      for i = 1:k
        temp = 1/( dr[i,di] + dl[k+1-i, di] )
        dtemp = Nderiv[i]*temp
        d2temp = N2deriv[i]*temp
        temp = bval[i]*temp

        N[i] = saved + dr[i,di]*temp
        Nderiv[i] = dsaved - temp + dr[i,di]*dtemp
        N2deriv[i] = d2saved - 2*dtemp + dr[i,di]*d2temp

        saved = dl[k+1-i,di]*temp
        dsaved = temp + dl[k+1-i,di]*dtemp
        d2saved = 2*dtemp + dl[k+1-i,di]*d2temp
      end  # End for i = 1:k
      N[k+1] = saved
      Nderiv[k+1] = dsaved
      N2deriv[k+1] = d2saved
    end # End for k = 1:order-1
  end  # End if order > 1

  return nothing
end
