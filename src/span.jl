# span.jl

@doc """
### findSpan

Determines the knot span index, i

**Inputs**

*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `U` : Knot vector
*  `k` : B-spline order
*  `nctl` : Number of control points

**Outputs**

*  `span` : Knot span index

SOURCE: The NURBS book 2nd Edition, Algorithm A2.1

"""->

function findSpan(u::Tffd, U::AbstractArray{Tffd,1}, k::Int, nctl::Int) where Tffd

  if u >= U[nctl+1]
    return nctl  # Special case when u = last term of knot vector
  elseif u < U[k]
    return k  # When u lies at the starting point of the curve
  end

  low = k-1
  high = nctl+1
  span = div(low + high, 2)
  # Do a binary search
  while u < U[span] || u >= U[span+1]
    if u < U[span]
      high = span
    else
      low = span
    end  # End if
    span = div(low+high, 2)
  end  # End while

  return span
end  # End function findSpan
