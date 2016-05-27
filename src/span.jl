# span.jl

@doc """
### findSpan

Determines the knot span index, i

**Inputs**

*  `u` : coordinate value (u,v,w is the coordinate space in 3D)
*  `map` : Object of type mapping
*  `di`  : coordinate direction in which `u` exists

**Outputs**

*  `span` : Knot span index

SOURCE: The NURBS book 2nd Edition, Algorithm A2.1

"""->

function findSpan(u, map, di)

  k = map.order[di]
  U = view(map.knot, :, di)

  if u >= U[n+1]
    return n  # Special case when u = last term of knot vector
  elseif u < U[k]
    return k  # When u lies at the starting point of the curve
  end

  low = k-1
  high = n+1
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
