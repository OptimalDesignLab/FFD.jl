# boundingbox.jl
# contain all functions pertaining to the bounding box

type BoundingBox

  # Physical space
  ndim::Integer # 2D or 3D
  geom_coord::Array{AbstractFloat, 2}  # highest and lowst values of x, y, z
                                       # coordinates along every dimension
                                       # in the physical space
  offset::Array{AbstractFloat, 1} # Offset from the mesh for the bounding box dimensions
  box_coord::Array{AbstractFloat, 2}
  lengths::Array{AbstractFloat, 1} # dimensions of the bounding box

  # Parametric space
  function BoundingBox(dim, coord, spacing) # 3D

    # Allocate members
    geom_coord = Array(AbstractFloat, 2, ndim) # each row = x_min, x_max,  each column = ndim
    box_coord = Array(AbstractFloat, 2, ndim)  # same as above
    offset = Array(AbstractFloat, ndim) # offset for each dimenstion in the (x, y, z) space
    lengths = Array(AbstractFloat, ndim) # same as above

    # Populate members of BoundingBox
    ndim = dim
    geom_coord[:,:] = coord[:,:]
    offset[:] = spacing[:]

    # Get the boumding box coordinates and dimensions
    for i = 1:size(geom_coord,2)
      box_coord[1,i] = geom_coord[1,i] - offset[i]
      box_coord[2,i] = geom_coord[2,i] + offset[i]
      lengths[i] = norm(box_coord[:,i], 2)
    end

    new(ndim, geom_coord, offset, box_coord, lengths)

  end
end
