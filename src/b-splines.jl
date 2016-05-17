# B-Spline Basis function for Bezier Bernsetirn splines
@doc """
### basisFunctions

Computes the basis functions needed to compute B-Splines

**Inputs**

*  `i` : Knot index
*  `u` : Knot
*  `p` : Degree of B-spline basis function (Order p+1)
*  `U` : Knot Vector
*  `N` : Array of basis functions

**Outputs**

*  None

SOURCE: The NURBS book 2nd Edition, Algorithm A2.2

"""->

function basisFunctions(i, u, p, U, N)

  # Initialize variables
  left = Array(Float64, p)
  right = Array(Float64, p)
  saved = 0.0

  N[1] = 1.0
  for j = 2:p+1
    left[j] = u - U[i+1-j]
    right[j] = U[i+j] - u
    saved = 0.0

    for r = 1:j
      temp = N[r]/(right[r+1] + left[j-r])
      N[r] = saved + right[r+1]*temp
      saved = left[r]*temp
    end

    N[j] = saved
  end

  return nothing
end
