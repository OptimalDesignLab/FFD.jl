# boundingbox.jl
# contain all functions pertaining to the bounding box

@doc """
###PumiBoundingBox

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

type PumiBoundingBox{Tffd} <: AbstractBoundingBox

  # Physical space
  ndim::Integer # 2D or 3D
  origin::Array{Tffd, 1}
  unitVector::Array{Tffd, 2}
  geom_bounds::Array{Tffd, 2}  # highest and lowst values of x, y, z
                                       # coordinates along every dimension
                                       # in the physical space
  offset::Array{Float64, 1} # Offset from the mesh for the bounding box dimensions
  box_bounds::Array{Tffd, 2}
  lengths::Array{Tffd, 1} # dimensions of the bounding box

  # Parametric space
  function PumiBoundingBox(map::PumiMapping, mesh::AbstractMesh,
                           sbp::AbstractSBP, offset) # 3D

    # Check Inputs
    for i = 1:3
      @assert offset[i] >= 0 "The bounding box offsets need to be >= 0 in every direction"
    end
    if mesh.dim == 2
      @assert offset[3] > 0 "2D meshes requires offset along 3rd dimension > 0"
    end

    ndim = 3 # The Bounding box is always 3D

    # Allocate members
    origin = zeros(Tffd, ndim)
    unitVector = eye(Tffd, ndim) # Create default a unit vector aligned with the physical space
    geom_bounds = zeros(Tffd, 2, ndim) # each row = x_min, x_max,  each column = ndim
    box_bounds = zeros(Tffd, 2, ndim)  # same as above
    lengths = zeros(Tffd, ndim) # same as above

    # Populate members of BoundingBox
    if map.full_geom == true
      if mesh.isDG == true
        calcEntireGeometryBounds(mesh.vert_coords, geom_bounds)
      else
        calcEntireGeometryBounds(mesh.coords, geom_bounds)
      end  # End if mesh.isDG == true
    else
      calcSurfaceGeomBounds(mesh, sbp, geom_bounds,map.geom_faces)
    end  # End if map.full_geom == true

    # Get the boumding box coordinates and dimensions
    for i = 1:size(geom_bounds,2)  # TODO: Come up with a better definition of box_bound and lengths
      box_bounds[1,i] = geom_bounds[1,i] - offset[i]
      box_bounds[2,i] = geom_bounds[2,i] + offset[i]
      lengths[i] = box_bounds[2,i] - box_bounds[1,i]
    end

    origin[:] = box_bounds[1,:] # Get the lower x,y,z ordinates to be defined as the
                               # origin

    return new(ndim, origin, unitVector, geom_bounds, offset, box_bounds, lengths)
  end
end

@doc """
### FreeFormDeformation.calcGeomBounds

Gets the geometric bounds of a pumi mesh object from its mesh coordinates. This
function works for both 2D and 3D pumi mesh objects

**Arguments**

* `coords` : Pumi Mesh coordinates
* `geom_bounds` : Highest and lowest values of x,y,z coordinates along every
                  dimension. dim1 = highest and lowest values, dim2 = x,y or z
                  direction i.e.
```
[xmin ymin zmin
 xmax ymax zmax]
```

"""->

function calcEntireGeometryBounds{Tffd}(coords::AbstractArray{Tffd,3},
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

  if !MPI.Initialized()
    MPI.Init()
  end
  recv_buffer = MPI.Allgather(geom_bounds, MPI.COMM_WORLD)
  for i = 1:6:length(recv_buffer)
    if recv_buffer[i] < geom_bounds[1,1]
      geom_bounds[1,1] = recv_buffer[i]
    end
    if recv_buffer[i+1] > geom_bounds[2,1]
      geom_bounds[2,1] = recv_buffer[i+1]
    end
    if recv_buffer[i+2] < geom_bounds[1,2]
      geom_bounds[1,2] = recv_buffer[i+2]
    end
    if recv_buffer[i+3] > geom_bounds[2,2]
      geom_bounds[2,2] = recv_buffer[i+3]
    end
    if recv_buffer[i+4] < geom_bounds[1,3]
      geom_bounds[1,3] = recv_buffer[i+4]
    end
    if recv_buffer[i+5] > geom_bounds[2,3]
      geom_bounds[2,3] = recv_buffer[i]
    end
  end

  return nothing
end

function calcSurfaceGeomBounds{Tffd}(mesh::AbstractCGMesh, sbp::AbstractSBP,
                                     geom_bounds::AbstractArray{Tffd,2},
                                     geom_faces::AbstractArray{Int,1})

  fill!(geom_bounds, 0.0)

  xmin = zero(Tffd)
  xmax = zero(Tffd)
  ymin = zero(Tffd)
  ymax = zero(Tffd)
  zmin = zero(Tffd)
  zmax = zero(Tffd)


  if mesh.dim == 2
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range)
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        for j = 1:sbp.numfacenodes
          k = sbp.facenodes[j, bndry_i.face]
          coords = view(mesh.coords, :, k, bndry_i.element)
          if xmin > coords[1]
            xmin = coords[1]
          elseif xmax < coords[1]
            xmax = coords[1]
          end # End if

          if ymin > coords[2]
            ymin = coords[2]
          elseif ymax < coords[2]
            ymax = coords[2]
          end # End if
        end   # End for j = 1:sbp.facenode
      end     # End for i = 1:nfaces
    end       # End for itr = 1:length(geom_faces)
  else
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end

      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:end_index
      bndry_facenums = view(mesh.bndryfaces, start_index:(end_index - 1))
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        for j = 1:sbp.numfacenodes
          k = sbp.facenodes[j, bndry_i.face]
          coords = view(mesh.coords, :, k, bndry_i.element)
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

          if zmin > coords[3,j,i]
            zmin = coords[3,j,i]
          elseif zmax < coords[3,j,i]
            zmax = coords[3,j,i]
          end # End if
        end   # End for j = 1:sbp.facenode
      end     # End for i = 1:nfaces
    end       # End for itr = 1:length(geom_faces)
  end

  geom_bounds[1,1] = xmin
  geom_bounds[2,1] = xmax
  geom_bounds[1,2] = ymin
  geom_bounds[2,2] = ymax
  geom_bounds[1,3] = zmin
  geom_bounds[2,3] = zmax

  if !MPI.Initialized()
    MPI.Init()
  end
  recv_buffer = MPI.Allgather(geom_bounds, MPI.COMM_WORLD)
  for i = 1:6:length(recv_buffer)
    if recv_buffer[i] < geom_bounds[1,1]
      geom_bounds[1,1] = recv_buffer[i]
    end
    if recv_buffer[i+1] > geom_bounds[2,1]
      geom_bounds[2,1] = recv_buffer[i+1]
    end
    if recv_buffer[i+2] < geom_bounds[1,2]
      geom_bounds[1,2] = recv_buffer[i+2]
    end
    if recv_buffer[i+3] > geom_bounds[2,2]
      geom_bounds[2,2] = recv_buffer[i+3]
    end
    if recv_buffer[i+4] < geom_bounds[1,3]
      geom_bounds[1,3] = recv_buffer[i+4]
    end
    if recv_buffer[i+5] > geom_bounds[2,3]
      geom_bounds[2,3] = recv_buffer[i]
    end
  end

  return nothing
end

function calcSurfaceGeomBounds{Tffd}(mesh::AbstractDGMesh, sbp::AbstractSBP,
                                     geom_bounds::AbstractArray{Tffd,2},
                                     geom_faces::AbstractArray{Int,1})

  fill!(geom_bounds, 0.0)

  xmin = zero(Tffd)
  xmax = zero(Tffd)
  ymin = zero(Tffd)
  ymax = zero(Tffd)
  zmin = zero(Tffd)
  zmax = zero(Tffd)


  if mesh.dim == 2
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          coords = view(mesh.vert_coords, :, vtx_arr[j], bndry_i.element)
          if xmin > coords[1]
            xmin = coords[1]
          elseif xmax < coords[1]
            xmax = coords[1]
          end # End if

          if ymin > coords[2]
            ymin = coords[2]
          elseif ymax < coords[2]
            ymax = coords[2]
          end # End if
        end   # End for j = 1:length(vtx_arr)
      end     # End for i = 1:nfaces
    end       # End for itr = 1:length(geom_faces)
  else
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],geom_face_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr2]
      end_index = mesh.bndry_offsets[itr2+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          coords = view(mesh.vert_coords, :, vtx_arr[j], bndry_i.element)
          if xmin > coords[1]
            xmin = coords[1]
          elseif xmax < coords[1]
            xmax = coords[1]
          end # End if

          if ymin > coords[2]
            ymin = coords[2]
          elseif ymax < coords[2]
            ymax = coords[2]
          end # End if

          if zmin > coords[3]
            zmin = coords[3]
          elseif zmax < coords[3]
            zmax = coords[3]
          end # End if
        end   # End for j = 1:length(vtx_arr)
      end     # End for i = 1:nfaces
    end       # End for itr = 1:length(geom_faces)
  end

  geom_bounds[1,1] = xmin
  geom_bounds[2,1] = xmax
  geom_bounds[1,2] = ymin
  geom_bounds[2,2] = ymax
  geom_bounds[1,3] = zmin
  geom_bounds[2,3] = zmax

  if !MPI.Initialized()
    MPI.Init()
  end
  recv_buffer = MPI.Allgather(geom_bounds, MPI.COMM_WORLD)
  for i = 1:6:length(recv_buffer)
    if recv_buffer[i] < geom_bounds[1,1]
      geom_bounds[1,1] = recv_buffer[i]
    end
    if recv_buffer[i+1] > geom_bounds[2,1]
      geom_bounds[2,1] = recv_buffer[i+1]
    end
    if recv_buffer[i+2] < geom_bounds[1,2]
      geom_bounds[1,2] = recv_buffer[i+2]
    end
    if recv_buffer[i+3] > geom_bounds[2,2]
      geom_bounds[2,2] = recv_buffer[i+3]
    end
    if recv_buffer[i+4] < geom_bounds[1,3]
      geom_bounds[1,3] = recv_buffer[i+4]
    end
    if recv_buffer[i+5] > geom_bounds[2,3]
      geom_bounds[2,3] = recv_buffer[i]
    end
  end

  return nothing
end
