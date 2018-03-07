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
#=
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
=#
#=
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
=#

"""
  Evalutes the new locations of the surface points.  Users should update
  `map.cp_xyz` and then call this function.

  **Inputs**

   * map: a PumiMapping object
   * mesh: a PumiMesh object

  **Inputs/Outputs**

   * pts: a mesh.dim x map.numFacePts array to be overwritten with the new
          surface points (ordered according to the numbering supplied to
          the `PumiMapping` constructor.
"""
function evalSurface{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractDGMesh, pts::Array{Tffd, 2})

  @assert size(pts, 1) == mesh.dim
  @assert size(pts, 2) == map.numFacePts

  x = zeros(Tffd, 3)  # compatability between 2D and 3D

  for i=1:map.numFacePts

    xi_i = sview(map.xi, :, i)
    evalVolumePoint(map, xi_i, x)

    for j=1:mesh.dim
      pts[j, i] = x[j]
    end
  end

  return nothing
end


#=
"""
  Computes the surface node coordinates based on the current control point
  positions.  Users should update the values in `map.cp_xyz` and then call this
  function.

  **Inputs**
  
   * map: the PumiMappping object
   * mesh: the Pumi mesh object

  **Outputs**

   * vertices: An array of arrays holding the new surface coordinates.
               The length of the outer array is the number of geometric
               entities the FFD box encapsulates.  Each inner array is
               `mesh.dim` x `mesh.coord_numNodesPerFace` x the number of faces
               on the geometric entity

  Note this function is currently limited to linear coordinate fields.
"""
function evalSurface{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractDGMesh)

  @assert mesh.coord_order == 1

  nwall_faces = getnWallFaces(mesh, map.bc_nums)
  vertices = Array(Array{Tffd,3}, length(map.bc_nums))
  defineVertices(mesh, map.bc_nums, vertices)

  # vertices = zeros(Tffd, size(mesh.vert_coords))
  x = zeros(Tffd, 3)
  for (idx, itr) in enumerate(map.bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      # get the local index of the vertices on the boundary face (local face number)
      # map.xi has the xi coordinates of the face nodes only
      for j = 1:mesh.coord_numNodesPerFace
        fill!(x, 0.0)
        evalVolumePoint(map, map.xi[idx][:,j,i], x)
        for k = 1:map.ndim
          vertices[idx][k, j, i] = x[k]
        end
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces
  end  # End for itr = 1:length(map.bc_nums)

  return vertices
end
=#
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

function evalVolumePoint{Tffd}(map::AbstractMappingType, xi::AbstractArray{Tffd,1}, xyz)

  fill!(xyz, 0.0)

  Nu = zeros(Tffd, map.order[1])
  Nv = zeros(Tffd, map.order[2])
  Nw = zeros(Tffd, map.order[3])

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

function contractWithdGdB{Tffd}(map::AbstractMappingType{Tffd}, xi, dJdG)

  # Evaluate the knot vectors and basis values
  Nu = zeros(Tffd, map.order[1])
  Nv = zeros(Tffd, map.order[2])
  Nw = zeros(Tffd, map.order[3])
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
        for idim = 1:3
          map.work[idim, startu+ii, startv+jj, startw+kk] += coeff*dJdG[idim]
        end
      end  # End for kk = 1:map.order[3]
    end    # End for jj = 1:map.order[2]
  end      # End for ii = 1:map.order[1]

  return nothing
end  # End function contractWithdGdB(map, dJdGrid)


"""
  Evalutes a transposed Jacobian-vector product

  Xcp_bar =  (dXs/dXcp)^T Xs_bar

  where Xs are the surface coordinates, and Xcp are the control point
  coordinates.

  **Inputs**

   * map: a PumiMapping object
   * mesh: a PumiMeshDG
   * Xs_bar: the array to multiply the transposed jacobian against,
             size 3 x `map.numFacePts`

  **Inputs/Outputs**

   * Xcp_bar: arrayto be overwritten with results, same size as `map.cp_xyz`
"""
function evaldXdControlPointTransposeProduct{Tmsh, T, Tffd}(map::PumiMapping{Tffd}, mesh::AbstractDGMesh{Tmsh}, Xs_bar::AbstractMatrix, Xcp_bar::AbstractArray{T, 4})

  @assert size(Xcp_bar, 1) == 3  # all meshes are 3 dimensional to FFD
  @assert size(Xcp_bar, 2) == map.nctl[1]
  @assert size(Xcp_bar, 3) == map.nctl[2]
  @assert size(Xcp_bar, 4) == map.nctl[3]

  @assert size(Xs_bar, 1) == mesh.dim
  @assert size(Xs_bar, 2) == map.numFacePts

  Xs_bar_i = zeros(Tffd, 3)  # compatability with 2D and 3D
  fill!(map.work, 0.0)
  for i=1:map.numFacePts
    xi_i = sview(map.xi, :, i)
    for j=1:mesh.dim
      Xs_bar_i[j] = Xs_bar[j, i]
    end
    contractWithdGdB(map, xi_i, Xs_bar_i)
  end

  work_slice = sview(map.work, 1:3, :, :, :)
  copy!(Xcp_bar, work_slice)

  return nothing
end

"""
  Evaluates the Jacobian-vector product

  Xs_dot = (dXs/dXcp) Xcp_dot

  where Xs are the surface coordinates, and Xcp are the control point
  coordinates.

  It is not safe to complex step this function.

  **Inputs**

   * map: a PumiMapping object
   * mesh: a PumiMeshDG
   * Xcp_dot: array to multiply the jacobian against, same size as `map.cp_xyz`

  **Inputs/Outputs**

   * Xs_dot: the array to be overwritten with the results
             size 3 x `map.numFacePts`

  **Implementation Notes**

  This function uses complex step internally to do the calculation, using
  `map.map_cs`.  `Tffd` does not have to be a complex type for the `map`
  passed to this function.


"""
function evaldXdControlPointProduct{Tmsh, T, Tffd}(map::PumiMapping{Tffd},
               mesh::AbstractDGMesh{Tmsh}, Xcp_dot::AbstractArray{T, 4},
               Xs_dot::AbstractMatrix )
#TODO: add some simd

  @assert size(Xs_dot, 1) == map.ndim
  @assert size(Xs_dot, 2) == map.numFacePts
  for i=1:4
    @assert size(Xcp_dot, i) == size(map.cp_xyz, i)
  end

  # copy things into map.map_cs
  map_cs = map.map_cs
  copy!(map_cs.cp_xyz, map.cp_xyz)
#  copy!(map_cs.xi, map.xi)

  h = 1e-20
  pert = Complex128(0, h)
  @simd for i=1:length(Xcp_dot)
    map_cs.cp_xyz[i] += Xcp_dot[i]*pert
  end

  vertices = map_cs.vertices
  evalSurface(map_cs, mesh, vertices)

  @simd for i=1:length(Xs_dot)
    Xs_dot[i] = imag(vertices[i])/h
  end

  # remove the perturbation
  @simd for i=1:length(Xcp_dot)
    map_cs.cp_xyz[i] -= Xcp_dot[i]*pert
  end

  return nothing
end



#=
"""
  evaluates (dXs/dXcp)^T v, where v is a user-supplied vector, Xs are the
  surface grid points and Xcp are the control points (note that Xs is called
  G and Xcp is called B in other parts of the code).

  Specifically, Xs is the non-unqiue representation of the surface vertices
  (it is per-geometric entity).

  **Inputs**

   * map: the PumiMapping object
   * mesh: the Pumi mesh object
   * Xs_bar: the v vector that the transposed jacobian is multiplied against.
             It is an array of arrays, where the outer array has length equal to
             the number of geometric entities the FFD box encapsulates.  Each
             inner array has dimensions `mesh.dim` x `mesh.coord_numNodesPerFace`
             x the number of faces on the geometric entity.

  **Inputs/Outputs**

   * Xcp_bar: array, same shape as map.cp_xyz.  This array will be overwritten
              with the result

  Note: Xs_bar allows duplicate entries (if a vertex is on the intersection of 
        two geometric entites), but only one entry can be non-zero
        in order to get the correct result
"""
function evaldXdControlPointTransposeProduct{Tmsh}(map::PumiMapping, mesh::AbstractDGMesh{Tmsh}, Xs_bar::Array{Array{Complex128, 3}, 1}, Xcp_bar::AbstractArray{Complex128, 4})

  @assert size(Xcp_bar, 1) == 3  # all meshes are 3 dimensional to FFD
  @assert size(Xcp_bar, 2) == map.nctl[1]
  @assert size(Xcp_bar, 3) == map.nctl[2]
  @assert size(Xcp_bar, 4) == map.nctl[3]

  @assert length(Xs_bar) == length(map.bc_nums)
  for i=1:length(Xs_bar)
    @assert size(Xs_bar[i], 1) == mesh.dim
    @assert size(Xs_bar[i], 2) == mesh.coord_numNodesPerFace
    @assert size(Xs_bar[i], 3) == size(map.xi[i], 3)
  end
  
  fill!(map.work, 0.0)  # prepare for accumulation
  for (idx, itr) in enumerate(map.bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)
    Xs_bar_itr = Xs_bar[idx]
    Xs_bar_j = zeros(Tmsh, 3)  # vector of length 3 gives compatability between
                               # 2D and 3D

    for i=1:nfaces
      for j=1:mesh.coord_numNodesPerFace
        # copy value into Xs_bar_j
        for k=1:mesh.dim
          Xs_bar_j[k] = Xs_bar_itr[k, j, i]
        end
        contractWithdGdB(map, map.xi[idx][:, j, i], Xs_bar_j)

        # contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
      end  # end loop j
    end  # end loop i
  end  # end loop itr

  # copy result into user-provided array
  # was 1:mesh.dim
  work_slice = sview(map.work, 1:3, :, :, :)
  copy!(Xcp_bar, work_slice)

  return nothing
end
=#








#=
"""
### evaldXdControlPointProduct

Higher level function that is used for computing the derivative

  ∂J/∂B = (∂J/∂G)(∂G/∂B)

for pumi meshes. Here J denotes the functional, G denotes the grid, and B denotes
the control point.

**Arguments**

* `map` : Object of PumiMapping type
* `mesh` : Object of AbstractDGMesh type
* `dJdVert` : derivative of the a functional w.r.t. the pumi mesh vertices.

"""
function evaldXdControlPointProduct(map::PumiMapping, mesh::AbstractDGMesh,
                                    dJdVert::AbstractArray{Float64,1})

  my_rank = MPI.Comm_rank(MPI.COMM_WORLD)

  fill!(map.work, 0.0)
  ctr = 1
  local_vertnum_history = Int[]
  dJdVert_arr = reshape(dJdVert, 3, convert(Int,length(dJdVert)/3))

  # Gather nface_verts from all ranks
  nlocal_face_verts = getLocalNumFaceVerts_unique(mesh, map.geom_faces)
  nlocal_face_verts_arr = MPI.Allgather(nlocal_face_verts, MPI.COMM_WORLD)
  ranks_for_bndry_faces = Int[]
  for i = 1:length(nlocal_face_verts_arr)
    if nlocal_face_verts_arr[i] > 0
      push!(ranks_for_bndry_faces, i-1)
    end # End if nface_verts_arr[i] > 0
  end

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
          # Check if this local vertex exists on other ranks
          if haskey(mesh.vert_sharing.rev_mapping, local_vertnum)
            # Check which all MPI ranks have this vertex and its position in
            # that rank's mesh.vert_sharing.
            rank_arr = first(mesh.vert_sharing.rev_mapping[local_vertnum])
            localIdx_arr = last(mesh.vert_sharing.rev_mapping[local_vertnum])
            # Lower rank always has ownership rights
            # if my_rank < minimum(rank_arr)
            #   contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
            #   ctr += 1
            # else
            #   mpi_neighbor_rank = getNeighborRank(mesh, bndry_i)
            #   if my_rank < mpi_neighbor_rank
            #     contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
            #     ctr += 1
            #   end
            # end # End if my_rank < minimum(rank_arr)
            intersect_arr = intersect(rank_arr, ranks_for_bndry_faces)
            if intersect_arr == []
              contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
              ctr += 1
            elseif my_rank < minimum(intersect_arr)
              contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
              ctr += 1
            end # End if intersect_arr == []
          else # The vertex exist only on one MPI rank
            contractWithdGdB(map, map.xi[itr][:,j,i], dJdVert_arr[:,ctr])
            ctr += 1
          end # End if haskey
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
=#
