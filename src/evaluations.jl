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

function evalVolume{Tffd}(map::PumiMapping{Tffd}, mesh::AbstractMesh)

  fill!(mesh.coords, 0.0)

  if mesh.dim == 2 # If its a 2D PumiMesh
    arr = zeros(Tffd, 3)
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        fill!(arr,0.0)
        evalVolumePoint(map, map.xi[:,j,i], arr)
        mesh.coords[:,j,i] = arr[1:2]
      end
    end
  else  # 3D pumi mesh
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        xyz = view(mesh.coords, :,j,i)
        evalVolumePoint(map, map.xi[:,j,i], xyz)
      end
    end
  end  # End If

  return nothing
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
    arr = zeros(Tffd, 3)
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
      itr2 = 0
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:end_index
      bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        for j = 1:sbp.numfacenodes
          fill!(x, 0.0)
          k = sbp.facenodes[j, bndry_i.face]
          arr[1:2] = mesh.coords[:,k,bndry_i.elements]
          evalVolumePoint(map, map.xi[:,j,i,itr], arr)
          mesh.coords[:,k,bndry_i.elements] = arr[1:2]
        end  # End for j = 1:sbp.numfacenodes
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
  else
    for itr = 1:length(map.geom_faces)
      geom_face_number = map.geom_faces[itr]
      itr2 = 0
      # get the boundary array associated with the geometric edge
      itr2 = 0
      for itr2 = 1:mesh.numBC
        if findfirst(mesh.bndry_geo_nums[itr2],g_edge_number) > 0
          break
        end
      end
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:end_index
      bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        for j = 1:sbp.numfacenodes
          fill!(x, 0.0)
          k = sbp.facenodes[j, bndry_i.face]
          xyz = view(mesh.coords,:,k,bndry_i.elements)
          evalVolumePoint(map, map.xi[:,j,i,itr], arr)
        end  # End for j = 1:sbp.numfacenodes
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
  end    # End if map.ndim == 2

  return nothing
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

function evalVolumePoint(map, xi, xyz)

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
  #=
  println("startu = $startu, startv = $startv, startw = $startw")
  println("Nu = \n$Nu")
  println("Nv = \n$Nv")
  println("Nw = \n$Nw")
  =#
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
