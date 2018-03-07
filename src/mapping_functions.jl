@doc """
### linearMap

Take in one coordinate in (x,y,z) space and convert it to parametric (s,t,u)
space. The symbols used in this function have been taken from the paper
"Free Form Deformation of Solid Geometric Models, Sederberg & Parry, 1986". It
may change at a later date

**Arguments**

*  `map` : Object of mapping type
*  `box` : BoundingBox object
*  `X`   : Point coordinate in (x,y,z) space
*  `pX`   : corresponding coordinates in (s,t,u) space

"""->

function linearMap(map::AbstractMappingType, box::AbstractBoundingBox,
                   X, pX)

  # The assumption of the mapping for the bounding box presently is that
  # coordinate transformation from physical x,y,z coordinate transformation to
  # parametric s,t,u coordinate involved only translation and scaling. There is
  # no rotation.

  # Get the x,y,z coordinates for the origin of the s,t,u system
  origin = box.origin
  S = box.unitVector[:,1]*box.lengths[1]
  T = box.unitVector[:,2]*box.lengths[2]
  U = box.unitVector[:,3]*box.lengths[3]

  XmX0 = X - origin

  # calculate s
  # Division by lengths[i] to ensure s,t,u lie between [0,1]
  TcrossU = cross(T,U)
  s = dot(TcrossU,XmX0)/dot(TcrossU,S)

  # Calculate t
  ScrossU = cross(S,U)
  t = dot(ScrossU,XmX0)/dot(ScrossU,T)

  # calculate u
  ScrossT = cross(S,T)
  u = dot(ScrossT,XmX0)/dot(ScrossT,U)

  pX[:] = [s,t,u]

  return nothing
end

@doc """
### nonlinearMap

Computes the (s,t,u) parametric coordinates of the a point in the geometry
embedded in the FFD box of any aribitrary shape. This is done by doing a Newton
solve.

**Arguments**

*  `map` : Object of Mapping type
*  `box` : bounding box object. In this case the bounding box can be any shape
*  `X`   : (x,y,z) coordinates of the point in geometry
*  `pX`  : (s,t,u) coordinates of the point. An initial guess of pX must be
           supplied. This is needed by the Newton's solve

"""->

function nonlinearMap{Tffd}(map::AbstractMappingType{Tffd},
                            box::AbstractBoundingBox{Tffd}, X, pX)

  origin = box.origin

  # Compute the residual
  res = zeros(Tffd, box.ndim)
  pointVal = zeros(Tffd, box.ndim)
  xi = zeros(Tffd, box.ndim)
  xi_new = zeros(xi)
  xi[:] = pX

  # Do the newton solve to get the (s,t,u coordinates)
  for itr = 1:50
    # Compute residual
    fill!(pointVal, 0.0)
    evalVolumePoint(map, xi, pointVal)
    res = X - pointVal
    # Construct jacobian
    J = zeros(box.ndim, box.ndim)
    jderiv = zeros(Int, box.ndim)

    for i = 1:box.ndim
      fill!(jderiv, 0)
      jderiv[i] = 1
      Jrow = view(J,i,:)
      calcdXdxi(map, xi, jderiv, Jrow)
    end
    xi_new = xi + J\res
    if norm(xi_new - xi, 2) < 1e-15
      pX[:] = xi_new[:]
      break
    else
      xi[:] = xi_new[:]
    end

  end

  return nothing
end

@doc """
### calcParametricMappingLinear

Creates a linear mapping for an array of nodes in the (x,y,z) space to the
(s,t,u) space.

**Arguments**

*  `map` : Object of Mapping type
*  `box` : BoundingBox object
*  `nodes_xyz` : (x,y,z) coordinates of the nodes of the mesh
"""->

function calcParametricMappingLinear(map::Mapping, box,
                                     nodes_xyz::AbstractArray{AbstractFloat,4})

  X = zeros(map.ndim)
  for k = 1:map.numnodes[3]
    for j = 1:map.numnodes[2]
      for i = 1:map.numnodes[1]
        X[:] = nodes_xyz[i,j,k,:]
        pX = view(map.xi,i,j,k,:)
        linearMap(map, box, X, pX)
      end
    end
  end

  return nothing
end  # End function calcParametricLinear

function calcParametricMappingLinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractCGMesh)

  @assert false "CG meshes not supported"

  if mesh.dim == 2
    X = zeros(Tffd,3)
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        X[1:2] = mesh.coords[:,j,i]
        pX = view(map.xi,:,j,i)
        linearMap(map, box, X, pX)
      end
    end
  else
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        X = view(mesh.coords,:,j,i)
        pX = view(map.xi,:,j,i)
        linearMap(map, box, X, pX)
      end
    end
  end

  return nothing
end

function calcParametricMappingLinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractDGMesh)

  @assert false "this function is broken"
  if mesh.dim == 2
    X = zeros(Tffd,3)
    for i = 1:mesh.numEl
      for j = 1:size(mesh.vert_coords,2) # 1:mesh.numNodesPerElement
        X[1:2] = mesh.vert_coords[:,j,i]
        pX = view(map.xi,:,j,i)
        linearMap(map, box, X, pX)
      end
    end
  else
    for i = 1:mesh.numEl
      for j = 1:size(mesh.vert_coords,2) # 1:mesh.numNodesPerElement
        X = view(mesh.vert_coords,:,j,i)
        pX = view(map.xi,:,j,i)
        linearMap(map, box, X, pX)
      end
    end
  end

  return nothing
end

function calcParametricMappingLinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractCGMesh,
                                     bc_nums::AbstractArray{Int,1})

  @assert false "CG meshes not supported"
  if mesh.dim == 2
    x = zeros(Tffd,3)
    for (idx, itr) in enumerate(bc_nums)
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:end_index
      bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          fill!(x, 0.0)
          x[1:2] = mesh.coords[:,vtx_arr[j],bndry_i.element]
          pX = sview(map.xi[idx], :, j, i)
          linearMap(map, box, x, pX)
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)
  else
    for (idx, itr) in enumerate(bc_nums)
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:end_index
      bndry_facenums = sview(mesh.bndryfaces, start_index:(end_index - 1))
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          fill!(x, 0.0)
          X = view(mesh.coords,:,vtx_arr,bndry_i.elements)
          pX = view(map.xi[idx], :, j, i)
          linearMap(map, box, X, pX)
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)

  end  # End if mesh.dim == 2

  return nothing
end

function calcParametricMappingLinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractDGMesh,
                                     bc_nums::AbstractArray{Int,1})

  @assert false "This function is broken"
  # Check if the knot vectors are for Bezier Curves with Bernstein polynomial
  # basis functions
  for i = 1:length(map.edge_knot)
    ctr = 0
    for j = 1:length(map.edge_knot[i])
      if map.edge_knot[i][j] != 0.0 && map.edge_knot[i][j] != 1.0
        ctr += 1
      end
    end
    @assert ctr == 0 "Linear mapping works only for Bezier Curves with Bernstein polynomial basis functions"
  end

  x = zeros(Tffd,3)
  for (idx, itr) in enumerate(bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      # get the local index of the vertices
      vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(vtx_arr)
        fill!(x, 0.0)
        for k = 1:mesh.dim
          x[k] = mesh.vert_coords[k,vtx_arr[j],bndry_i.element]
        end
        pX = view(map.xi[idx], :, j, i)
        linearMap(map, box, x, pX)
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(geomfaces)

  return nothing
end

@doc """
### calcParametricMappingNonlinear

Creates a non linear mapping for an array of nodes in the (x,y,z) space to the
(s,t,u) space.

**Arguments**

*  `map` : Object of Mapping type
*  `box` : BoundingBox object
*  `nodes_xyz` : (x,y,z) coordinates of the points in the embedded geometry

"""->

function calcParametricMappingNonlinear(map::Mapping, box,
                                        nodes_xyz::AbstractArray{AbstractFloat,4})

  X = zeros(map.ndim)
  pX = zeros(map.ndim)
  for k = 1:map.numnodes[3]
    for j = 1:map.numnodes[2]
      for i = 1:map.numnodes[1]
        X[:] = nodes_xyz[i,j,k,:]
        pX[:] = [0.,0.,0.]
        nonlinearMap(map, box, X, pX)
        map.xi[i,j,k,:] = pX[:]
      end
    end
  end

  return nothing
end

# unused now that full_geom is no more?
function calcParametricMappingNonlinear{Tffd}(map::PumiMapping{Tffd},
                                        box::PumiBoundingBox, mesh::AbstractCGMesh)

  @assert false "CG meshes not supported"
  if mesh.dim == 2
    X = zeros(Tffd,3)
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        X[1:2] = mesh.coords[:,j,i]
        pX = view(map.xi,:,j,i)
        nonlinearMap(map, box, X, pX)
      end
    end
  else
    for i = 1:mesh.numEl
      for j = 1:mesh.numNodesPerElement
        X = view(mesh.coords,:,j,i)
        pX = view(map.xi,:,j,i)
        nonlinearMap(map, box, X, pX)
      end
    end
  end

  return nothing
end

# unused now that full_geom is no more?
function calcParametricMappingNonlinear{Tffd}(map::PumiMapping{Tffd},
                                        box::PumiBoundingBox, mesh::AbstractDGMesh)

  @assert false "This function is broken"
  X = zeros(Tffd,3)
  for i = 1:mesh.numEl
    for j = 1:mesh.numNodesPerElement
      for k = 1:mesh.dim
        X[k] = mesh.vert_coords[k,j,i]
      end
      pX = view(map.xi,:,j,i)
      nonlinearMap(map, box, X, pX)
    end
  end

  return nothing
end

function calcParametricMappingNonlinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractCGMesh,
                                     bc_nums::AbstractArray{Int,1})
  @assert false "CG Meshes not supported"
  if mesh.dim == 2
    x = zeros(Tffd,3)
    for (idx, itr) in enumerate(bc_nums)
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range)
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          fill!(x, 0.0)
          x[1:2] = mesh.coords[:,vtx_arr[j],bndry_i.element]
          pX = view(map.xi[idx], :, j, i)
          nonlinearMap(map, box, x, pX)
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)
  else
    for (idx, itr) in enumerate(bc_nums)
      start_index = mesh.bndry_offsets[itr]
      end_index = mesh.bndry_offsets[itr+1]
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range)
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        for j = 1:length(vtx_arr)
          X = view(mesh.coords,:,vtx_arr[j],bndry_i.element)
          pX = view(map.xi, :, j, i)  
          nonlinearMap(map, box, X, pX)
        end  # End for j = 1:length(vtx_arr)
      end    # End for i = 1:nfaces
    end      # End for itr = 1:length(geomfaces)

  end  # End if mesh.dim == 2

  return nothing
end

function calcParametricMappingNonlinear{Tffd}(map::PumiMapping{Tffd},
                                     box::PumiBoundingBox, mesh::AbstractDGMesh,
                                     bc_nums::AbstractArray{Int,1})

  x_real = zeros(Float64, 3)
  for i=1:map.numFacePts
    v_i = map.face_verts[i]

    getPoint(mesh.m_ptr, v_i, 0, x_real)
    pX = sview(map.xi, :, i)
    nonlinearMap(map, box, sview(x_real, :), pX)
  end
#=
  for (idx, itr) in enumerate(bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range) # faces on geometric edge i
    nfaces = length(bndry_facenums)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      # get the local index of the vertices
      vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
      for j = 1:length(vtx_arr)
        fill!(x, 0.0)
        for k = 1:mesh.dim
          x[k] = mesh.vert_coords[k,vtx_arr[j],bndry_i.element]
        end
        pX = view(map.xi[idx], :, j, i)
        nonlinearMap(map, box, x, pX)
      end  # End for j = 1:length(vtx_arr)
    end    # End for i = 1:nfaces
  end      # End for itr = 1:length(geomfaces)

  return nothing
end
=#
  return nothing
end
@doc """
### calcdXdxi

It calculates the partial derivative of the mapping (including the 0th order)

(x,y,z) = [x(xi,eta,zeta), y(xi,eta,zeta), z(xi,eta,zeta)]

**Arguments**

*  `map`    : Object of Mapping type
*  `xi`     : Mapping parametric coordinates (s,t,u) to be evaluated
*  `jderiv` : Derivative indices. jderi[mdi] = n means take the nth derivative
              in the xi[mdi] direction.
*  `dX`     : Derivative of X w.r.t the 3 parametric variables. length = 3

REFERENCE: Carl de Boor, 'Pratical Guide to Splines', pgs 138-149, function
           BVALUE

Notes: (taken from Mapping_Mod.f90)
1) de Boor points out (pg 149) that for many points (as we may
   have) it is faster to compute the piecewise polys and differentiate
   them.  Something to consider for the future.
2) Since direction indices are used for both the physical and
   mathematical coordinates, di will be reserved for the physical
   and mdi for the mathematical coordinates in this function.

"""->

function calcdXdxi(map, xi, jderiv, dX)

  # initialize intermediate arrays
  km1 = zeros(Int, 3)
  jcmax = zeros(Int, 3)
  jcmin = zeros(Int, 3)
  imk = zeros(Int, 3)
  left = zeros(Int, 3)

  # the derivative is zero if any of jderiv(mdi) >= order(mdi)
  for mdi = 1:3
    if jderiv[mdi] >= map.order[mdi]
      return
    end
  end

  # find the spatially varying knot vector
  # calcKnot(map)

  # find the left(3) array such that
  #   knot(mdi,left(mdi)) <= xi(mdi) <= knot(mdi,left(mdi)+1)
  for mdi = 1:3
    left[mdi] = findSpan(xi[mdi], map.edge_knot[mdi], map.order[mdi], map.nctl[mdi])
  end
  # the calculations involving the knot vectors are independent,
  # so loop through them sequentially


    for mdi = 1:3

      # we will store the (order) b-spline coefficients relevant
      # to the knot interval [knot(mdi,left),knot(mdi,left+1)]
      # in aj[mdi,1],...,aj[mdi,order]
      # and compute dl[mdi,j] = xi[mdi] - knot[mdi][left-j]
      #             dr[mdi,j] = knot[mdi][left+j] - xi[mdi]
      # for all j = 1,...,order[mdi]-1

      km1[mdi] = map.order[mdi] - 1
      jcmin[mdi] = 1
      imk[mdi] = left[mdi] - map.order[mdi]
      if imk[mdi] < 0
        # we are close to the left end of the knot interval, so some
        # of the aj will be set to zero later
        jcmin[mdi] = 1 - imk[mdi]
        for j = 1:left[mdi]
          map.dl[j,mdi] = xi[mdi] - map.edge_knot[mdi][left[mdi]+1-j]
        end
        for j = left[mdi]:km1[mdi]
          map.dl[j,mdi] = map.dl[left[mdi],mdi]
        end
      else
        for j = 1:km1[mdi]
          map.dl[j,mdi] = xi[mdi] - map.edge_knot[mdi][left[mdi]+1-j]
        end
      end

      jcmax[mdi] = map.order[mdi]
      n = map.nctl[mdi]
      nmi = n - left[mdi]
      if nmi < 0
        # we are close to the right end of the knot interval, so some
        # of the aj will be set to zero later
        jcmax[mdi] = map.order[mdi] + nmi
        for j = 1:jcmax[mdi]
          map.dr[j,mdi] = map.edge_knot[mdi][left[mdi]+j] - xi[mdi]
        end
        for j = jcmax[mdi]:km[mdi]
          map.dr[j,mdi] = map.dr[jcmax[mdi],mdi]
        end
      else
        for j = 1:km1[mdi]
          map.dr[j,mdi] = map.edge_knot[mdi][left[mdi]+j] - xi[mdi]
        end
      end
    end # mdi loop
    # set all elements of aj[:,:,1] to zero, in case we are close to
    # the ends of the knot vector
    map.aj[:,:,1] = 0.0

    for jc1 = jcmin[1]:jcmax[1]
      p = imk[1] + jc1

      # set all the elements of aj[:,:,2] to zero, in case we are
      # close to the ends of the knot vector
      map.aj[:,:,2] = 0.0

      for jc2 = jcmin[2]:jcmax[2]
        q = imk[2] + jc2

        # set all the elements of aj[:,:,3] to zero, in case we are
        # close to the ends of the knot vector
        map.aj[:,:,3] = 0.0

        for jc3 = jcmin[3]:jcmax[3]
          map.aj[:,jc3,3] = map.cp_xyz[:,p,q,imk[3]+jc3]
        end

        if jderiv[3] != 0
          # derivative: apply the recursive formula X.12b from de Boor
          for j = 1:jderiv[3]
            kmj = map.order[3] - j
            ilo = kmj
            for jj = 1:kmj
              map.aj[:,jj,3] = (map.aj[:,jj+1,3] - map.aj[:,jj,3]) *
                               kmj / (map.dl[ilo,3] + map.dr[jj,3])
              ilo = ilo - 1
            end
          end
        end
        #println("map.aj[:,:,3] = \n", map.aj[:,:,3])
        if jderiv[3] != km1[3]
          # if jderiv != order - 1, we need to apply the recursive
          # formula from de Boor
          for j = jderiv[3]+1:km1[3]
            kmj = map.order[3] - j
            ilo = kmj
            for jj = 1:kmj
              map.aj[:,jj,3] = (map.aj[:,jj+1,3]*map.dl[ilo,3] +
                   map.aj[:,jj,3]*map.dr[jj,3]) / (map.dl[ilo,3] + map.dr[jj,3])
              ilo = ilo - 1
            end
          end
        end
        # println("map.aj[:,:,3] = \n", map.aj[:,:,3])

        map.aj[:,jc2,2] = map.aj[:,1,3]
      end # get next element of aj(:,2)

      if jderiv[2] != 0
        # derivative: apply the recursive formula X.12b from de Boor
        for j = 1:jderiv[2]
          kmj = map.order[2] - j
          ilo = kmj
          for jj = 1:kmj
            map.aj[:,jj,2] = (map.aj[:,jj+1,2] - map.aj[:,jj,2]) * kmj /
                             (map.dl[ilo,2] + map.dr[jj,2])
            ilo = ilo - 1
          end
        end
      end

      if jderiv[2] != km1[2]
        # if jderiv /= order - 1, we need to apply the recursive
        # formula from de Boor
        for j = jderiv[2]+1:km1[2]
          kmj = map.order[2] - j
          ilo = kmj
          for jj = 1:kmj
            map.aj[:,jj,2] = (map.aj[:,jj+1,2]*map.dl[ilo,2] + map.aj[:,jj,2]*
                             map.dr[jj,2]) / (map.dl[ilo,2] + map.dr[jj,2])
            ilo = ilo - 1
          end
        end
      end

      map.aj[:,jc1,1] = map.aj[:,1,2]
      # println("map.aj[:,:,2] = \n", map.aj[:,:,2])

    end # get next element of map.aj[:,:,1]

    if jderiv[1] != 0
      # derivative: apply the recursive formula X.12b from de Boor
      for j = 1:jderiv[1]
        kmj = map.order[1] - j
        ilo = kmj
        for jj = 1:kmj
          map.aj[:,jj,1] = (map.aj[:,jj+1,1] - map.aj[:,jj,1]) * kmj /
                           (map.dl[ilo,1] + map.dr[jj,1])
          ilo = ilo - 1
        end
      end
    end

    if jderiv[1] != km1[1]
      # if jderiv /= order - 1, we need to apply the recursive
      # formula from de Boor
      for j = jderiv[1]+1:km1[1]
        kmj = map.order[1] - j
        ilo = kmj
        for jj = 1:kmj
          map.aj[:,jj,1] = (map.aj[:,jj+1,1]*map.dl[ilo,1] + map.aj[:,jj,1]*
                           map.dr[jj,1]) / (map.dl[ilo,1] + map.dr[jj,1])
          ilo = ilo - 1
        end
      end
    end
    # println("map.aj[:,:,1] = \n", map.aj[:,:,2])

    dX[:] = map.aj[:,1,1]

  return nothing
end  # End function calcdXdxi(map, xi, jderiv)
