# boundingbox.jl
# contain all functions pertaining to the bounding box

@doc """
### BoundingBox

Class that contains information necessary to create a bonding box of control
points.

**Members**

*  `ndim` : 2D or 3D
*  `origin` : intended origin of the (s,t,u) parametric coordinate system in the
              (x,y,z) space
*  `unitVector`: Unit vectors that define the (s,t,u) coordinate axes expressed
                 in (x,y,z) coordinate system. dim1 = x,y,z values, dim2 = s,t
                 or u vector
*  `geom_coord` : Highest and lowest values of x,y,z coordinates along every
                  dimension. dim1 = highest and lowest values, dim2 = x,y or z
                  direction
*  `offset` : Offset from the geometry/mesh boundary for determining the size of
              the bounding box
*  `box_bound` : Ordinates of the extremas of the bounding box. dim1 = upper
                 and lower bounds, dim2 = x,y or z direction.
*  `lengths` : lengths of the edges of the parallelopied being generated

The current assumption behind this type is that the parametric axes has the same
orientation as the physical coordinate system.

"""->

type BoundingBox

  # Physical space
  ndim::Integer # 2D or 3D
  origin::Array{AbstractFloat, 1}
  unitVector::Array{AbstractFloat, 2}
  geom_coord::Array{AbstractFloat, 2}  # highest and lowst values of x, y, z
                                       # coordinates along every dimension
                                       # in the physical space
  offset::Array{AbstractFloat, 1} # Offset from the mesh for the bounding box dimensions
  box_bound::Array{AbstractFloat, 2}
  lengths::Array{AbstractFloat, 1} # dimensions of the bounding box

  # Parametric space
  function BoundingBox(dim, coord, spacing) # 3D

    ndim = dim # Set dimensions

    # Allocate members
    origin = Array(AbstractFloat, 3)
    unitVector = eye(ndim) # Create a unit vector aligned with the physical space
    geom_coord = Array(AbstractFloat, 2, ndim) # each row = x_min, x_max,  each column = ndim
    box_bound = Array(AbstractFloat, 2, ndim)  # same as above
    offset = Array(AbstractFloat, ndim) # offset for each dimenstion in the (x, y, z) space
    lengths = Array(AbstractFloat, ndim) # same as above

    # Populate members of BoundingBox
    geom_coord[:,:] = coord[:,:]
    offset[:] = spacing[:]

    # Get the boumding box coordinates and dimensions
    for i = 1:size(geom_coord,2)  # TODO: Come up with a better definition of box_bound and lengths
      box_bound[1,i] = geom_coord[1,i] - offset[i]
      box_bound[2,i] = geom_coord[2,i] + offset[i]
      # println("box_bound[:,i] = ", box_bound[:,i])
      lengths[i] = box_bound[2,i] - box_bound[1,i]
    end

    origin = box_bound[1,:] # Get the lower x,y,z ordinates to be defined as the
                            # origin

    new(ndim, origin, unitVector, geom_coord, offset, box_bound, lengths)

  end
end
