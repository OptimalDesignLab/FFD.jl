# knot.jl
# File for creating knots in the physical space in 1D, 2D or 3D.
@doc """
### calcKnot

Computes a uniform knot vector for the mapping object. It is assumed that the
number of control points is known from before. It is user specified. The knot
vector is computed in the parametric space [0,1] in all the 3 directions

**Inputs**

*  `map` : Object of mapping type
*  `xi`  : Parametric coordinate values

**Outputs**

*  None

"""->

function calcKnot(map, xi)

  for di = 1:map.ndim
    order = map.order[di]
    nctl = map.nctl[di]
    # Zero out the knot vector
    map.edge_knot[di][:] = 0.0
    step = 1 / (nctl-order+1)
    low = 0.0
    for i = (order+1):nctl
      map.edge_knot[di][i] = low + step
      low += step
    end
    # Fill the end of the knot vector
    map.edge_knot[di][nctl+1:nctl+order] = 1.0
  end

  return nothing
end

@doc """
### calcKnotLMS

Computes the knot vector in the parametric space [0,1]. This is for a least
squares curve approximation.

REFERENCE: Equation 9.68 & eqation 9.69, Page 412 The NURBS book, 2nd Edition

**Inputs**

*  `order` : Order of B-spline
*  `nctl`  : Number of control points for that knot vector
*  `X`     : ordinates used to make the knot vector in parametric space
*  `U`     : Knot vector to be populated

**Outputs**

*  None

"""->

function calcKnotLMS(order, nctl, X, U)
  # REFERENCE : Gaetan's pySpline knots.f90, knots_lms

  @assert length(U) = nctl + order "Length of specified knot vector is not equal
  to the sum of number of control points and order"

  lengthX = length(X)

  # Populate the extremas of the knot vector
  for i = 1:order
    U[i] = 0.0
    U[nctl+i] = 1.0
  end  # End for i = 1:order

  # Populate the middle of the knot vector. Recall the structure of knot vectors
  if isodd(lengthX)
    d = lengthX/(nctl-order+1)
    for j = 1:nctl-order
      i = convert(Int, floor(j*d)) # Conversion for integer based indexing
      alpha = j*d-i
      U[order+j] = (1-alpha)*X[i] + alpha*X[i+2]
    end
  else
    d = lengthX/(nctl-order+1)
    for j = 1:nctl-order
      i = convert(Int, floor(j*d)) # Conversion for integer based indexing
      alpha = j*d-i+0.5
      U[order+j] = (1-alpha)*X[i] + alpha*X[i]
    end
  end  # End if isodd(lengthX)

  return nothing
end

@doc """
### calcKnotInterp

Create a knot vector for an interpolating spline. This is again done in the
parametric space between [0,1]

**Inputs**

*  `X`  : Ordinates used to make the knot vector in the parametric space
*  `nd` : Number of derivatives specified
*  `order` : Order of the B-spline curve
*  `U`  : Knot vector to be computed

**Outputs**

*  None

"""->
function calcKnotInterp(X, nd, order, U)

  lengthX = length(X)
  nctl = lengthX + nd

  @assert nctl+order == length(U) # Sanity check

  # Populate the ends of the knot vector in the parametric space
  for j = 1:order
    U[j] = 0.0
    U[nctl+j] = 1.0
  end

  if nd == lengthX
    if order == 3
      for i = 1:(lengthX-2)
        U[order+2*i-1] = 0.5*(X[i] + X[i+1])
        U[order+2*i] = X[i+1]
      end # end for i = 1:(lengthX-2)
      U[nctl] = 0.5*(X[lengthX-1] + 1)
    elseif k == 4
      U[5] = 0.5 * X[2]
      U[nctl] = 0.5 * (X[lengthX-1] + 1)
      for i = 1:(lengthX-3)
        U[order+2*i] = (1/3)*(2*X[i+1] + X[i+2])
        U[order+2*i+1] = (1/3)*(X[i+1] + 2*X[i+2])
      end  # End for i = 1:(lengthX-3)
    end  # end if order == 3
  elseif nd == 0
    if order == 2
      for j = 1:lengthX-order
        U[order+j] = X[j+1]
      end
    elseif order == 3
      for j = 1:lengthX-order
        U[order+j] = 0.5*(X[j+1] + X[j+2])
      end
    elseif order == 4
      for j = 1:lengthX-order
        U[order+j] = X[j+2]
      end
    else
      error("Interpolation is only available for k = 2, 3 or 4")
    end
  else
    error("More options for nd not available yet")
  end  # End if nd == lengthX

  return nothing
end
