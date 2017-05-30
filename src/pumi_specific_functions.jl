# test_functions.jl
# File containing functions that reduce the suze of the test files

@doc """
###FreeFormDeformation.getnWallFaces

Get the number of element faces that exist on a geometric face of a mesh
despite the name of the function

**Inputs**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed

**Output**

* `nwall_faces` : Number of element faces(edges in 2D) on each geometric face
"""->

function getnWallFaces(mesh::AbstractDGMesh, geom_faces::AbstractArray{Int,1})

  # Prepare the wall coordinates array for mesh warping
  nwall_faces = zeros(Int,length(geom_faces))
  vtx_per_face = mesh.dim # only true for simplex elements
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
    nwall_faces[itr] = length(bndry_facenums)
  end      # End for itr = 1:length(geomfaces)

  return nwall_faces
end

@doc """
###FreeFormDeformation.getUniqueVertexArray

Get a 2D array of coordinates of unique vertices for the LOCAL portion of the
entire mesh.

**Input**

* `mesh` : Pumi DG mesh

**Output**

* `volNodes` : 2D array of unique vertices that exist on an MPI rank.
               size = [3, n_local_vertices]

"""->

function getUniqueVertexArray{Tmsh}(mesh::AbstractMesh{Tmsh})

  volNodes = zeros(Tmsh, 3, mesh.numVert)
  for i = 1:mesh.numEl
    for j = 1:size(mesh.vert_coords,2)
      # Get the vertex numbering on the portion of mesh owned by the processor
      local_vertnum = mesh.element_vertnums[j,i]
      for k = 1:mesh.dim
        volNodes[k, local_vertnum] = mesh.vert_coords[k,j,i] # mesh.element_vertnums
      end
    end
  end

  return volNodes
end

@doc """
###FreeFormDeformation.getUniqueWallCoordsArray

Get a 2D array of coordinates of unique vertices for the local portion of a
geometric face.

NOTE: Use this ONLY to compute the unique wall vertices on a given MPI rank.
      This function does not account for non-uniqueness that arises from the
      same vertex existing on 2 MPI ranks

**Inputs**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed

**Outputs**

* `wallCoords` : 2D array of coordinates of wall vertices. size = [3, total_vertices]

"""->

function getUniqueWallCoordsArray{Tmsh}(mesh::AbstractMesh{Tmsh},
                                  geom_faces::AbstractArray{Int,1})

  vtx_arr = getUniqueVertexArray(mesh) # Get Unique vertex Array
  wallCoords = getUniqueWallCoordsArray(mesh, vtx_arr, geom_faces)

  return wallCoords
end

function getUniqueWallCoordsArray{Tmsh}(mesh::AbstractMesh{Tmsh},
                                  vtx_arr::AbstractArray{Tmsh,2},
                                  geom_faces::AbstractArray{Int,1})

  nface_verts = getLocalNumFaceVerts_unique(mesh, geom_faces)
  wallCoords = zeros(Tmsh, 3, nface_verts)

  local_vertnum_history = Int[]
  ctr = 1
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
      # get the local index of the vertices
      face_vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(face_vtx_arr)
        local_vertnum = mesh.element_vertnums[face_vtx_arr[j],bndry_i.element]
        if findfirst(local_vertnum_history, local_vertnum) == 0
          for k = 1:mesh.dim
            wallCoords[k,ctr] = vtx_arr[k, local_vertnum]
          end
          ctr += 1
        end
        push!(local_vertnum_history, local_vertnum)
      end     # for j = 1:length(face_vtx_arr)
    end       # for i = 1:nfaces
  end # End for itr = 1:length(geom_faces)

  return wallCoords
end

@doc """
###FreeFormDeformation.getGlobalUniqueWallCorrdsArray

Get a 2D array of coordinates of unique wall vertices that is owned by the MPI
rank. By default PUMI Mesh does not have ownership considering the unstructured 
nature of the mesh where a vertex exists on several ranks. A sense of ownership
is created where the lowest rank containing the vertex owns it. The Array is 
still distributed across the MPI ranks however the shared vertices appear on
only the lowest rank.

**Input**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed

**Output**

* `wallCoords` : 2D array of wall coordinates owned by the rank. 
                 size = [3, n_owned_vertices]

"""->

function getGlobalUniqueWallCorrdsArray{Tmsh}(mesh::AbstractMesh{Tmsh},
                                  geom_faces::AbstractArray{Int,1})

  comm = MPI.COMM_WORLD
  comm_world = MPI.MPI_COMM_WORLD
  comm_self = MPI.COMM_SELF
  my_rank = MPI.Comm_rank(comm)
  comm_size = MPI.Comm_size(comm)

  wallCoords = zeros(Tmsh, 3, 0)
  vertex_coordinate = zeros(Tmsh,3)
  vtx_arr = getUniqueVertexArray(mesh) # Get Unique vertex Array

  nface_verts = getLocalNumFaceVerts_unique(mesh, geom_faces)
  if nface_verts != 0
    local_vertnum_history = Int[]
    for itr = 1:length(geom_faces)
      geom_face_number = geom_faces[itr]
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
        face_vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(face_vtx_arr)
          local_vertnum = mesh.element_vertnums[face_vtx_arr[j],bndry_i.element]
          if findfirst(local_vertnum_history, local_vertnum) == 0
            # Now Check if this local vertex exists on multiple ranks
            if haskey(mesh.vert_sharing.rev_mapping, local_vertnum)

              # Check which all MPI ranks have this vertex and its position in
              # that rank's mesh.vert_sharing.
              rank_arr = first(mesh.vert_sharing.rev_mapping[local_vertnum])
              localIdx_arr = last(mesh.vert_sharing.rev_mapping[local_vertnum])
              # The Idea is that the lowest rank always has ownership rights
              if my_rank < minimum(rank_arr)
                for k = 1:mesh.dim
                  vertex_coordinate[k] = vtx_arr[k, local_vertnum]
                end # End for k = 1:mesh.dim
                # Append a colum of coordinates
                wallCoords = hcat(wallCoords, vertex_coordinate)
              end # End if my_rank == minimum(rank_arr)

            else # The vertex exist only on one MPI rank

              for k = 1:mesh.dim
                vertex_coordinate[k] = vtx_arr[k, local_vertnum]
              end
              # Append a colum of coordinates
              wallCoords = hcat(wallCoords, vertex_coordinate)

            end   # End if haskey(mesh.vert_sharing.rev_mapping, local_vertnum)
          end     # End if findfirst(local_vertnum_history, local_vertnum) == 0
          push!(local_vertnum_history, local_vertnum)
        end     # for j = 1:length(face_vtx_arr)
      end       # for i = 1:nfaces
    end # End for itr = 1:length(geom_faces)
  end # End if nface_verts != 0

  MPI.Barrier(comm)

  return wallCoords
end

@doc """
###FreeFormDeformation.getLocalNumFaceVerts_unique

Obtain the number of unique face vertices that exist on one MPI rank. This
function DOES NOT take into account non-uniquness of vertices that exist across
multiple MPI ranks.

**Inputs**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed

**Output** :

* `ctr` : Number of unique vertices on a particular MPI rank
"""->
function getLocalNumFaceVerts_unique(mesh::AbstractDGMesh, geom_faces::AbstractArray{Int,1})

  local_vertnum_history = Int[]
  ctr = 0
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
      # get the local index of the vertices
      face_vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(face_vtx_arr)
        local_vertnum = mesh.element_vertnums[face_vtx_arr[j],bndry_i.element]
        if findfirst(local_vertnum_history, local_vertnum) == 0
          ctr += 1
        end
        push!(local_vertnum_history, local_vertnum)
      end # End for j = 1:length(face_vtx_arr)
    end   # End for i = 1:nfaces
  end

  return ctr
end

@doc """
###FreeFormDeformation.defineVertices

This function defines the shape of the vertices array that is used to update
the Pumi mesh after FFD when a geometric face is paramtereized.

**Arguments**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed
* `vertices` : Array of arrays holding the coodinates of the updated vertices
               shape = vertices[n_geom_faces][mesh.dim, size(mesh.vert_coords,2), n_elem_faces]

"""->

function defineVertices(mesh::AbstractDGMesh, geom_faces::AbstractArray{Int,1},
                        vertices::AbstractArray)

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
    vertices[itr] = zeros(size(mesh.vert_coords,1), size(mesh.vert_coords,2), nfaces)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      vertices[itr][:,:,i] = mesh.vert_coords[:,:,bndry_i.element]
    end # End for i = 1:nfaces
  end

  return nothing
end

@doc """
###FreeFormDeformation.getWallCoords

Get a 2D array of coordinates that exist on a wall element face. This function
does not take into account uniqueness across adjacent element faces that share
a common vertex or across MPI ranks that share a vertex.

**Inputs**

* `mesh` : Pumi DG mesh
* `geom_faces` : Array of geometric faces (edges in 2D) over which the number 
                 of element faces needs to be computed

**Output**

* `wallCoords` : 2D array of coordinates of vertices that exist on a wall.
                 size = [3, sum(nwall_faces)*vtx_per_face]

"""->

function getWallCoords(mesh::PumiMeshDG2, geom_faces::AbstractArray{Int, 1})

  nwall_faces = getnWallFaces(mesh, geom_faces)

  # Populate wallCoords
  vtx_per_face = mesh.dim # only true for simplex elements
  nWallCoords = sum(nwall_faces)*vtx_per_face
  wallCoords = zeros(3, nWallCoords)
  ctr = 1 # Counter for wall coordinates
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
      # get the local index of the vertices
      vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(vtx_arr)
        wallCoords[1:2, ctr] = mesh.vert_coords[:,vtx_arr[j],bndry_i.element]
        ctr += 1
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces
  end

  return wallCoords
end # End function getWallCoords
