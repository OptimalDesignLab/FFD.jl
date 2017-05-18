# test_functions.jl
# File containing functions that reduce the suze of the test files

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

function getUniqueWallCoordsArray{Tmsh}(mesh::AbstractMesh{Tmsh},
                                  geom_faces::AbstractArray{Int,1},
                                  is_face_closed::Bool=false)

  vtx_arr = getUniqueVertexArray(mesh) # Get Unique vertex Array
  wallCoords = getUniqueWallCoordsArray(mesh, vtx_arr, geom_faces, is_face_closed)

  return wallCoords
end

function getUniqueWallCoordsArray{Tmsh}(mesh::AbstractMesh{Tmsh},
                                  vtx_arr::AbstractArray{Tmsh,2},
                                  geom_faces::AbstractArray{Int,1},
                                  is_face_closed::Bool)

  nwall_faces = getnWallFaces(mesh, geom_faces)
  if is_face_closed == false
    wallCoords = zeros(Tmsh, 3, sum(nwall_faces) + 1)
  else
    wallCoords = zeros(Tmsh, 3, sum(nwall_faces))
  end

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

function checkValueExistence(Arr, val)

  ctr = 0
  for i = 1:length(Arr)
    if Arr[i] == val
      ctr += 1
    end
  end

  if ctr > 0
    return true
  else
    return false
  end

end
#=
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
=#
