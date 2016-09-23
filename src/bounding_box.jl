# boundingbox.jl
# contain all functions pertaining to the bounding box

@doc """
### BoundingBox

Class that contains information necessary to create a bonding box of control
points. Use the inner constructor to create an object of this type. The
arguments to the constructor in sequence are as follows

* Number of dimensions
* Highest and lowest values of x,y,z coordinates along every dimension.
  dim1 = highest and lowest values, dim2 = x,y or z direction
* Offset from the geometry/mesh boundary for determining the size of the
  bounding box

The Definition of the function is subject to change in the future. Presently,
the unit vectors that define the parametric coordinates of the box are assumed
to be aligned with the physical coordinate system. They can be changed
separately if needed. There is also an implicit assumption that the box is an
othogonal parallelopiped, but, this assumption is not enforce anywhere. The
final shape of the BoundingBox object will depend on the geometric location of
the control points in the physical space.

**Members**

*  `ndim` : 2D or 3D
*  `origin` : intended origin of the (s,t,u) parametric coordinate system in the
              (x,y,z) space
*  `unitVector`: Unit vectors that define the (s,t,u) coordinate axes expressed
                 in (x,y,z) coordinate system. dim1 = x,y,z values, dim2 = s,t
                 or u vector
*  `geom_bound` : Highest and lowest values of x,y,z coordinates along every
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

type BoundingBox{Tffd} <: AbstractBoundingBox

  # Physical space
  ndim::Integer # 2D or 3D
  origin::Array{Tffd, 1}
  unitVector::Array{Tffd, 2}
  geom_bound::Array{Tffd, 2}  # highest and lowst values of x, y, z
                                       # coordinates along every dimension
                                       # in the physical space
  offset::Array{Tffd, 1} # Offset from the mesh for the bounding box dimensions
  box_bound::Array{Tffd, 2}
  lengths::Array{Tffd, 1} # dimensions of the bounding box

  # Parametric space
  function BoundingBox(dim, coord, spacing) # 3D

    ndim = dim # Set dimensions

    # Allocate members
    origin = Array(AbstractFloat, 3)
    unitVector = eye(ndim) # Create default a unit vector aligned with the physical space
    geom_bound = Array(AbstractFloat, 2, ndim) # each row = x_min, x_max,  each column = ndim
    box_bound = Array(AbstractFloat, 2, ndim)  # same as above
    offset = Array(AbstractFloat, ndim) # offset for each dimenstion in the (x, y, z) space
    lengths = Array(AbstractFloat, ndim) # same as above

    # Populate members of BoundingBox
    geom_bound[:,:] = coord[:,:]
    offset[:] = spacing[:]

    # Get the boumding box coordinates and dimensions
    for i = 1:size(geom_bound,2)  # TODO: Come up with a better definition of box_bound and lengths
      box_bound[1,i] = geom_bound[1,i] - offset[i]
      box_bound[2,i] = geom_bound[2,i] + offset[i]
      lengths[i] = box_bound[2,i] - box_bound[1,i]
    end

    origin[:] = box_bound[1,:] # Get the lower x,y,z ordinates to be defined as the
                            # origin

    new(ndim, origin, unitVector, geom_bound, offset, box_bound, lengths)

  end
end

@doc """


"""->

function calcGeomBounds{Tffd}(coords::AbstractArray{Tffd,3},
                              geom_bounds::AbstractArray{Tffd,2})

  fill!(geom_bounds, 0.0)

  xmin = zero(Tffd)
  xmax = zero(Tffd)
  ymin = zero(Tffd)
  ymax = zero(Tffd)
  zmin = zero(Tffd)
  zmax = zero(Tffd)

  if size(coords,1) == 2
    for i = 1:size(coords,3)
      for j = 1:size(coords,2)
        if xmin > coords[1,j,i]
          xmin = coords[1,j,i]
        elseif xmax < coords[1,j,i]
          xmax = coords[1,j,i]
        end # End if

        if ymin > coords[2,j,i]
          ymin = coords[2,j,i]
        elseif ymax < coords[2,j,i]
          ymax = coords[2,j,i]
        end # End if
      end   # End for j = 1:size(coords,2)
    end     # End for i = 1:size(coords,3)
  elseif size(coords,1) == 3
    for i = 1:size(coords,3)
      for j = 1:size(coords,2)
        if xmin > coords[1,j,i]
          xmin = coords[1,j,i]
        elseif xmax < coords[1,j,i]
          xmax = coords[1,j,i]
        end

        if ymin > coords[2,j,i]
          ymin = coords[2,j,i]
        elseif ymax < coords[2,j,i]
          ymax = coords[2,j,i]
        end

        if zmin > coords[3,j,i]
          zmin = coords[3,j,i]
        elseif zmax < coords[3,j,i]
          zmax = coords[3,j,i]
        end # End if
      end   # End for j = 1:size(coords,2)
    end     # End for i = 1:size(coords,3)
  end

  geom_bounds[1,1] = xmin
  geom_bounds[2,1] = xmax
  geom_bounds[1,2] = ymin
  geom_bounds[2,2] = ymax
  geom_bounds[1,3] = zmin
  geom_bounds[2,3] = zmax

  return nothing
end
