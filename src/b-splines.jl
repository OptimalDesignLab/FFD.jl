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
function basisFunctions(map, N, di, u, span)

  # (U, order, u, span, N)
  order = map.order[di]
  U = view(map.knot, :, di)

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

Computes the basis functions and its first & second (method 2) derivative at a
particular point u.

**Inputs**

*  `map` : Object of mapping type
*  `N`   : Basis functions
*  `Nderiv` : First derivative of the basis functions at u
*  `N2deriv`: Second derivative of the basis functions at u
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
  dl = view(map.dl, :, :)

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

function derivBasisFunctions(map, N, Nderiv, N2deriv, di, u, span)

  order = map.order[di]
  N[1] = 1.0
  Nderiv[1] = 0.0
  N2deriv[1] = 0.0

  dr = view(map.dr, :, :)
  dl = view(map.dl, :, :)

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
      end  E End for i = 1:k
      N[k+1] = saved
      Nderiv[k+1] = dsaved
      N2deriv[k+1] = d2saved
    end # End for k = 1:order-1
  end  # End if order > 1

  return nothing
end

#=
@doc """
### basisFunctions3D

It returns the non-zero basis values for a 3D B-spline tensor product. This a 3D
version of basisFunctions.

**Arguments**

*  `map`    : Object of mapping type
*  `span`   : Array of Knot span indices. length = 3
*  `uvw`    : Parameter values
*  `bval3D` : 3D basis functions. each dimension corresponding to the dimensions
              of uvw. length = max_order

REFERENCE : Carl de Boor, 'Practical Guide to Splines', pgs 132-134

"""->

function basisFunctions3D(map, span, uvw, bval3D)

  # Use a Cartesian tensor product to find the basis values

  # Starting with di = 3 and work down
  order = map.order[3]
  offset = span[3]
  if order > 1
    for k = 1:(order-1)
      map.dr[k,3] = map.knot[offset+k,3] - uvw[3]
      map.dl[k,3] = uvw[3] - map.knot[offset+1-k,3]
      saved = 0.0
      for i = 1:k
        temp = bval3D[:,:,i] ./ (map.dr[i,3] + map.dl[k+1-i,3])
        bval3D[:,:,i] = saved + map.dr[i,3]*temp
        saved = map.dl[k+1-i,3]*temp
      end  # End for i = 1:k
      bval3D[:,:,k+1] = saved
    end  # end for k = 1:(order-1)
  end # end if order > 1

  # di = 2
  order = map.order[2]
  offset = span[2]
  if order > 1
    for k = 1:(order-1)
      map.dr[k,2] = map.knot[offset+k,2] - uvw[2]
      map.dl[k,2] = uvw[2] - map.knot[offset+1-k,2]
      saved = 0.0
      for i = 1:k
        temp = bval3D[:,i,:] ./ (map.dr[i,2] + map.dl[k+1-i,2])
        bval3D[:,i,:] = saved + map.dr[i,2]*temp
        saved = map.dl[k+1-i,2]*temp
      end  # End for i = 1:k
      bval3D[:,k+1,:] = saved
    end  # end for k = 1:(order-1)
  end # end if order > 1

  # di = 1
  order = map.order[1]
  offset = span[1]
  if order > 1
    for k = 1:(order-1)
      map.dr[k,1] = map.knot[offset+k,1] - uvw[1]
      map.dl[k,1] = uvw[1] - map.knot[offset+1-k,1]
      saved = 0.0
      for i = 1:k
        temp = bval3D[i,:,:] ./ (map.dr[i,1] + map.dl[k+1-i,1])
        bval3D[i,:,:] = saved + map.dr[i,1]*temp
        saved = map.dl[k+1-i,1]*temp
      end  # End for i = 1:k
      bval3D[k+1,:,:] = saved
    end  # end for k = 1:(order-1)
  end # end if order > 1

  @assert (sum(bval3D) - 1) < 1e-13

  return nothing
end
=#
