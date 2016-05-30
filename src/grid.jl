# Grid.jl

@doc """
### calcGrid

Calculates the mesh coordinates based on the element node values in the map,
hermite basis functions and the refinement level. The basis function values of
map should be known

**Arguments**

*  `map` : Object of Mapping type
*  `xyz` : Array of mesh node coordinates in j,k,m,1:3 format

"""->

function calcGrid(map, xyz)

  @assert size(map.jkmmax) == size(xyz)

  for m = 1:map.jkmmax[3]
    for k = 1:map.jkmmax[2]
      for j = 1:map.jkmmax[1]
        xi = view(map.ki, j, k, m, :)
        dXdxi(map, xi, jderiv, xyz[j,k,m,:])
      end  # End for j = 1:map.jkmmax[1]
    end  # End for k = 1:map.jkmmax[2]
  end  # End for m = 1:map.jkmmax[3]

  return nothing
end
