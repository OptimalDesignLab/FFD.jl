# Grid.jl

####  -------------------

#  NOT NEEDED !!!!!!!!!!

#### ---------------------

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


@doc """
### fitGrid

Finds the control point values such that, when evaluated at in computational
space, the mapping produces the mesh xyz as close as possible. Parameter
correction was added by L. Olague and J. Hicken and ported to the language
`Julia`

**Arguments**

*  `map` : Object of Mapping type
*  `xyz` : Mesh to be approximated. On exit, htis is over written with the new
           mesh produced by the map
*  `L2_error` : The L2 error squared in fitting the grid is added to this
                variable at the end of routine
*  `max_error`: The values in this array are updated if the infinity norm in the
                fit of this grid is larger

"""


function fitGrid(map, xyz, L2_error, max_error)

  # Allocate and initialize arrays


  return nothing
end
