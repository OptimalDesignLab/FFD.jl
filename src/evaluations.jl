# Evaluations

@doc """
### evalCurve

Determine the value of curve at certain points

**Inputs**

*  `u` : Array of coordinates where the curve needs to be computed
*  `U` : Knot vector
*  `order` : order of B-spline basis functions, (order = p+1)
*  `P` : Array of control points
*  `C` : Resulting curve values

**Outputs**

*  None

SOURCE: The NURBS book 2nd Edition, Algorithm A3.1
      : Gaetan's pyspline/src/eval_curve
"""->

function evalCurve(u, U, order, P, C)

  @assert length(P) + order == length(U)
  p = order - 1  # Degree of B-spline basis function
  nctl = length(P)
  N = Array(Float64, order) # Array of basis functions 1D

  for i = 1:length(u)
    span = findSpan(u[i], U, order, nctl)
    basisFunctions(U, order, u[i], span, N)
    C[i] = 0.0
    for j = 1:order
      C[i] += N[j]*P[span-order+j]
    end  # End for j = 1:p+1
  end    # End for i = 1:length(u)

  return nothing
end

@doc """
### evalVolume

Determine the the coordinates of all points in a 3D volume using B-splines. The
symbol convention used in the function is from the book

"The NURBS book 2nd Edition"

**Arguments**

*  `map` : Object of mapping type
*  `Vol` : (x,y,z) coordinates of the embedded volume within the contol points

####Method 2

Applies only to DG meshes so fae

**Input**

* `map`
* `mesh` : Abstract Pumi DG Mesh

**Output**

* `vertices` : Array of updated vertices same shape as mesh.vert_coords

"""->

function evalVolume{Tmsh}(map::Mapping, Vol::AbstractArray{Tmsh,4})

  fill!(Vol, 0.0) # Zero out all entries of Vol

  for k = 1:map.numnodes[3]
    for j = 1:map.numnodes[2]
      for i = 1:map.numnodes[1]
        xyz = view(Vol, i,j,k,:)
        evalVolumePoint(map, map.xi[i,j,k,:], xyz)
      end  # End for i = 1:map.numnodes[3]
    end    # End for j = 1:map.numnodes[2]
  end      # End for k = 1:map.numnodes[1]

  return nothing
end

function evalVolume{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractDGMesh)

  vertices = zeros(Tffd, size(mesh.vert_coords))
  arr = zeros(Tffd, 3)
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      fill!(arr,0.0)
      evalVolumePoint(map, map.xi[:,j,i], arr)
      for k = 1:map.ndim
        vertices[k,j,i] = arr[k]
      end
    end
  end

  return vertices
end


@doc """
evalSurface

Evaluate surface points (edge points for 2D) of a a pumi mesh object

**Arguments**

*  `map`  : PumiMapping object
*  `mesh` : Pumi mesh object
*  `sbp`  : Summation-By-Parts object

"""->

function evalSurface{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractCGMesh,
                           sbp::AbstractSBP)

  if map.ndim == 2
    x = zeros(Tffd, 3)
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
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
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          fill!(x, 0.0)
          evalVolumePoint(map, map.xi[itr][:,j,i], x)
          mesh.coords[:,vtx_arr[j],bndry_i.element] = x[1:2]
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
  else
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
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
      bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          xyz = view(mesh.coords,:,vtx_arr[j],bndry_i.element)
          evalVolumePoint(map, map.xi[itr][:,j,i], arr)
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
  end    # End if map.ndim == 2

  return nothing
end

function evalSurface{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractDGMesh)

  nwall_faces = getnWallFaces(mesh, map.geom_faces)
  vertices = Array(Array{Tffd,3}, length(map.geom_faces))
  defineVertices(mesh, map.geom_faces, vertices)

  # vertices = zeros(Tffd, size(mesh.vert_coords))
  x = zeros(Tffd, 3)
  for itr = 1:length(map.geom_faces)
    geom_face_number = map.geom_faces[itr]
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
        fill!(x, 0.0)
        evalVolumePoint(map, map.xi[itr][:,j,i], x)
        for k = 1:map.ndim
          vertices[itr][k,vtx_arr[j],i] = x[k]
        end
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces
  end  # End for itr = 1:length(map.geom_faces)

  return vertices
end

@doc """
### evalVolumePoint

Computes a point in the FFD volume. The symbol convention used in this function
is from "The NURBS book 2nd Edition".

**Arguments**

*  `map` : Object of Mapping type
*  `xi`  : Parametric FFD coordinates (referred to as the (s,t,u) coordinate
           systems within FFD functions and as (u,v,w) in this function)
*  `xyz` : (x,y,z) coordinates of the parametric point in the FFD volume
"""->

function evalVolumePoint(map::AbstractMappingType, xi, xyz)

  fill!(xyz, 0.0)

  Nu = zeros(map.order[1])
  Nv = zeros(map.order[2])
  Nw = zeros(map.order[3])

  # Work with u
  span = findSpan(xi[1], map.edge_knot[1], map.order[1], map.nctl[1])
  basisFunctions(map.edge_knot[1], map.order[1], xi[1], span, Nu)
  startu = span - map.order[1]

  # Work with v
  span = findSpan(xi[2], map.edge_knot[2], map.order[2], map.nctl[2])
  basisFunctions(map.edge_knot[2], map.order[2], xi[2], span, Nv)
  startv = span - map.order[2]

  # Work with w
  span = findSpan(xi[3], map.edge_knot[3], map.order[3], map.nctl[3])
  basisFunctions(map.edge_knot[3], map.order[3], xi[3], span, Nw)
  startw = span - map.order[3]

  for ii = 1:map.order[1]
    for jj = 1:map.order[2]
      for kk = 1:map.order[3]
        for idim = 1:3
          xyz[idim] += Nu[ii]*Nv[jj]*Nw[kk]*
                       map.cp_xyz[idim, startu+ii, startv+jj, startw+kk]
        end
      end  # End for kk = 1:map.order[3]
    end    # End for jj = 1:map.order[2]
  end      # End for ii = 1:map.order[1]

  return nothing
end

@doc """
### contractWithdGdB

Used to contract the matrix dG/dB, the derivative of the mesh node coordinates
with respect to the mapping control points, with a given array, here labelled
dJdGrid and stored in (j,k,m,:) format. This can be used to find
dJ/dB = (dJ/dG)(dG/dB) for some objective J. In practice, the objective used is
actually

             J_practice = J_obj + psi^T dot R

where J_obj is the actual objective, psi are the flow adjoints, and R is the
flow residual.

**Arguments**

*  `map`     : Object of mapping type
*  `dJdGrid` : The derivative of the objective w.r.t the grid coordintes

"""->

function contractWithdGdB(map::AbstractMappingType, xi, dJdG)

  # Evaluate the knot vectors and basis values
  Nu = zeros(map.order[1])
  Nv = zeros(map.order[2])
  Nw = zeros(map.order[3])
  span = zeros(Int, 3)

  # Work with u
  span[1] = findSpan(xi[1], map.edge_knot[1], map.order[1], map.nctl[1])
  basisFunctions(map.edge_knot[1], map.order[1], xi[1], span[1], Nu)
  startu = span[1] - map.order[1]

  # Work with v
  span[2] = findSpan(xi[2], map.edge_knot[2], map.order[2], map.nctl[2])
  basisFunctions(map.edge_knot[2], map.order[2], xi[2], span[2], Nv)
  startv = span[2] - map.order[2]

  # Work with w
  span[3] = findSpan(xi[3], map.edge_knot[3], map.order[3], map.nctl[3])
  basisFunctions(map.edge_knot[3], map.order[3], xi[3], span[3], Nw)
  startw = span[3] - map.order[3]

  for ii = 1:map.order[1]
    for jj = 1:map.order[2]
      for kk = 1:map.order[3]
        coeff = Nu[ii]*Nv[jj]*Nw[kk]
        # println("ii = $ii, jj = $jj, kk = $kk, startu = $startu, startv = $startv, startw = $startw")
        for idim = 1:3
          map.work[idim, startu+ii, startv+jj, startw+kk] += coeff*dJdG[idim]
        end
      end  # End for kk = 1:map.order[3]
    end    # End for jj = 1:map.order[2]
  end      # End for ii = 1:map.order[1]

  return nothing
end  # End function contractWithdGdB(map, dJdGrid)

function evaldXdControlPointProduct(map::PumiMapping, mesh::AbstractDGMesh,
                                    dJdVert::AbstractArray{Float64,1})

  fill!(map.work, 0.0)
  ctr = 1
  local_vertnum_history = Int[]
  dJdVert_arr = reshape(dJdVert, 3, convert(Int,length(dJdVert)/3))
  for itr = 1:length(map.geom_faces)
    geom_face_number = map.geom_faces[itr]
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
      vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(vtx_arr)
        local_vertnum = mesh.element_vertnums[vtx_arr[j],bndry_i.element]
        if findfirst(local_vertnum_history, local_vertnum) == 0
          contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
          ctr += 1
        end
        push!(local_vertnum_history, local_vertnum)
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces

    # Now do an All reduce operation on map.work
    recv_arr = zeros(map.work)
    MPI.Allreduce!(map.work, recv_arr, MPI.SUM, MPI.COMM_WORLD)
    for i = 1:length(map.work)
      map.work[i] = recv_arr[i]
    end

  end  # End for itr = 1:length(map.geom_faces)

  return nothing
end
