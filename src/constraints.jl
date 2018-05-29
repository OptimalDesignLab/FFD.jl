# This file contains functions for setting optimization constraints on the FFD
# mapping.

@doc """
### numLinearPlaneConstraints

Returns the number of plane constraints (see setLinearPlaneConstraints! for
further details).

**Arguments**

* `map` : PumiMapping object
* `di` : the index the defines the constrained plane (e.g. `di`=2 -> j planes)

"""->
function numLinearPlaneConstraints(map::PumiMapping{Tffd}, di::Int) where Tffd
  it1 = mod(di,3)+1
  it2 = mod(di+1,3)+1
  num_per_plane = (map.nctl[it1]*map.nctl[it2] - 3)*3
  return map.nctl[di]*num_per_plane
end

@doc """
### countVarsLinearPlaneConstraints!

Updates array `cntEq` with the number of variables involved in each linear-plane
constraint.

**Arguments**

* `map` : PumiMapping object
* `di` : the index the defines the constrained plane (e.g. `di`=2 -> j planes)
* `cntEq` : array with number of variables involved in each equality constraint
* `fncidx` : index to the next equality constraint to add to cntEq

**Returns** `fncidx` + number of linear-plane constraints + 1

"""->
function countVarsLinearPlaneConstraints!(map::PumiMapping{Tffd}, di::Int,
                                          cntEq::AbstractArray{Int,1},
                                          fncidx::Int) where Tffd
  num_cnstr = numLinearPlaneConstraints(map, di)
  cntEq[fncidx:fncidx+num_cnstr-1] = 4
  return fncidx+num_cnstr
end

@doc """
### setLinearPlaneConstraints!

Set constraints that relate CPs in given planes to CPs at the corners of the
plane.  Basically, a slice of CPs in the FFD are constrained to remain in a
plane, and shear and translate according to the CPs at the corners of the plane.

**Arguments**

* `map`  : PumiMapping object
* `di`   : the index the defines the constrained plane (e.g. `di`=2 -> j planes)
* `iLfun`, `jLvar` : the k^th element in the gradient matrix has row index
*                    `iLfun(k)` and column index `jLvar(k)`
* `LinG` : `LinG(k) = A(iLfun(k),jLvar(k))`
* `Flow`, `Fupp` : the upper and lower bounds on all constraints
* `fncidx` : function counter; next function to use in iLfun
* `ptr` : index pointer to next free space in iLfun and jLvar

**Returns** `fncidx` + number of linear-plane constraints + 1, and `ptr` + total
  number of non-zeros in the sparse contraints + 1

"""->
function setLinearPlaneConstraints!(map::PumiMapping{Tffd}, di::Int, iLfun::AbstractArray{Int,1},
    jLvar::AbstractArray{Int,1}, LinG::AbstractArray{Tffd,1},
    Flow::AbstractArray{Tffd,1}, Fupp::AbstractArray{Tffd,1},
    fncidx::Int, ptr::Int) where Tffd

  it1 = mod(di,3)+1
  it2 = mod(di+1,3)+1

  jkm = zeros(Int,3)
  cp1_jkm = zeros(jkm); cp2_jkm = zeros(jkm); cp3_jkm = zeros(jkm)
  cp1_idx = zeros(jkm); cp2_idx = zeros(jkm); cp3_idx = zeros(jkm)
  cp1_xyz = zeros(3); cp2_xyz = zeros(3); cp3_xyz = zeros(3)
  u1 = zeros(3); u2 = zeros(3); u3 = zeros(3); e1 = zeros(3); e2 = zeros(3)
  x = zeros(3)

  # loop over the CP planes in the direction di
  for jdi = 1:map.nctl[di]
    jkm[di] = jdi

    # set the (j,k,m) indices for the CPs at the corners of the plane
    cp1_jkm[di] = jdi
    cp1_jkm[it1] = 1
    cp1_jkm[it2] = 1

    cp2_jkm[di] = jdi
    cp2_jkm[it1] = 1
    cp2_jkm[it2] = map.nctl[it2]

    cp3_jkm[di] = jdi
    cp3_jkm[it1] = map.nctl[it1]
    cp3_jkm[it2] = 1

    # get the coordinates of the CPs at the corners of the plane, and some unit
    # vectors
    cp1_xyz = map.cp_xyz[:,cp1_jkm[1],cp1_jkm[2],cp1_jkm[3]]
    cp2_xyz = map.cp_xyz[:,cp2_jkm[1],cp2_jkm[2],cp2_jkm[3]]
    cp3_xyz = map.cp_xyz[:,cp3_jkm[1],cp3_jkm[2],cp3_jkm[3]]

    u1[:] = cp2_xyz - cp1_xyz
    fac1 = 1./norm(u1)
    u1[:] *= fac1

    u2[:] = cp3_xyz - cp1_xyz
    fac2 = 1./norm(u2)
    u2[:] *= fac2

    u3[:] = cross(u1,u2)
    u3[:] *= 1./norm(u3)
    e1[:] = cross(u2,u3)
    e2[:] = cross(u3,u1)

    # get the indices of the CPs at the corner of the plane
    cp1_idx = map.cp_idx[:,cp1_jkm[1],cp1_jkm[2],cp1_jkm[3]]
    cp2_idx = map.cp_idx[:,cp2_jkm[1],cp2_jkm[2],cp2_jkm[3]]
    cp3_idx = map.cp_idx[:,cp3_jkm[1],cp3_jkm[2],cp3_jkm[3]]

    # loop over the CPs that lie in the jdi-th plane
    for jit1 = 1:map.nctl[it1]
      jkm[it1] = jit1
      for jit2 = 1:map.nctl[it2]
        jkm[it2] = jit2
        if jkm == cp1_jkm || jkm == cp2_jkm || jkm == cp3_jkm
          # the corner CPs are independent, so skip them
          continue
        end
        # get the dependency coefficients
        x[:] = map.cp_xyz[:,jkm[1],jkm[2],jkm[3]] - cp1_xyz
        u1_coeff, u2_coeff = calcObliqueProjection(x, u1, u2, e1, e2)
        coeff2 = u1_coeff*fac1
        coeff3 = u2_coeff*fac2
        coeff1 = -coeff2 - coeff3

        # insert into sparse storage
        for bdi = 1:3
          # set data for dependent CP coord
          iLfun[ptr] = fncidx
          jLvar[ptr] = map.cp_idx[bdi,jkm[1],jkm[2],jkm[3]]
          LinG[ptr] = 1.0
          ptr += 1
          # set data for independent CP1
          iLfun[ptr] = fncidx
          jLvar[ptr] = cp1_idx[bdi]
          LinG[ptr] = -coeff1
          ptr += 1
          # set data for independent CP2
          iLfun[ptr] = fncidx
          jLvar[ptr] = cp2_idx[bdi]
          LinG[ptr] = -coeff2
          ptr += 1
          # set data for independent CP3
          iLfun[ptr] = fncidx
          jLvar[ptr] = cp3_idx[bdi]
          LinG[ptr] = -coeff3
          ptr += 1
          Flow[fncidx] = 0.0
          Fupp[fncidx] = 0.0
          fncidx += 1
        end

      end # jit2 loop
    end # jit1 loop
  end # jdi loop
  return fncidx, ptr
end

@doc """
### calcObliqueProjection

Basis coefficients for u1 and u2 based on an oblique projection along e1, e2.

**Arguments**

* `x` : vector whose representation in (`u1`,`u2`) is sought
* `u1`,`u2` : unit basis vectors, not necessarily orthogonal
* `e1`,`e2` : covariant basis vectors satisfying dot(e1,u2) = 0, dot(e2,u1) = 0

"""->
function calcObliqueProjection(x::AbstractArray{T,1}, u1::AbstractArray{T,1},
                               u2::AbstractArray{T,1}, e1::AbstractArray{T,1},
                               e2::AbstractArray{T,1}) where T
  return dot(e1,x)/dot(e1,u1), dot(e2,x)/dot(e2,u2)
end

@doc """
### numLinearCornerConstraints

Returns the number of corner constraints (see setLinearCornerConstraints! for
further details).

**Arguments**

* `map` : PumiMapping object
* `di` : the index the defines the constrained plane (e.g. `di`=2 -> j planes)

"""->
function numLinearCornerConstraints(map::PumiMapping{Tffd}, di::Int) where Tffd
  return map.nctl[di]*2
end

@doc """
### countVarsLinearCornerConstraints!

Updates array `cntEq` with the number of variables involved in each linear-corner
constraint.

**Arguments**

* `map` : PumiMapping object
* `di` : the index the defines the constrained plane (e.g. `di`=2 -> j planes)
* `cntEq` : array with number of variables involved in each equality constraint
* `fncidx` : index to the next equality constraint to add to cntEq

**Returns** `fncidx` + number of linear-plane constraints + 1

"""->
function countVarsLinearCornerConstraints!(map::PumiMapping{Tffd}, di::Int,
                                           cntEq::AbstractArray{Int,1},
                                           fncidx::Int) where Tffd
  num_cnstr = numLinearCornerConstraints(map, di)
  cntEq[fncidx:fncidx+num_cnstr-1] = 4
  return fncidx+num_cnstr
end

@doc """
### setLinearCornerConstraints!

This constraint is used to keep each plane in the `di` direction square; this
should be used in conjuction with setLinearPlaneConstraints!.  It also assumes
that the di coordinates of the CPs are the same.

**Arguments**

* `map`  : PumiMapping object
* `di`   : the index the defines the constrained plane (e.g. `di`=2 -> j planes)
* `iLfun`, `jLvar` : the k^th element in the gradient matrix has row index
*                    `iLfun(k)` and column index `jLvar(k)`
* `LinG` : `LinG(k) = A(iLfun(k),jLvar(k))`
* `Flow`, `Fupp` : the upper and lower bounds on all constraints
* `fncidx` : function counter; next function to use in iLfun
* `ptr` : index pointer to next free space in iLfun and jLvar

**Returns** `fncidx` + number of linear-plane constraints + 1, and `ptr` + total
  number of non-zeros in the sparse contraints + 1

"""->
function setLinearCornerConstraints!(map::PumiMapping{Tffd}, di::Int, iLfun::AbstractArray{Int,1},
    jLvar::AbstractArray{Int,1}, LinG::AbstractArray{Tffd,1},
    Flow::AbstractArray{Tffd,1}, Fupp::AbstractArray{Tffd,1},
    fncidx::Int, ptr::Int) where Tffd

  it1 = mod(di,3)+1
  it2 = mod(di+1,3)+1

  jkm = zeros(Int,3)
  cp1_jkm = zeros(jkm); cp2_jkm = zeros(jkm); cp3_jkm = zeros(jkm)
  cp1_idx = zeros(jkm); cp2_idx = zeros(jkm); cp3_idx = zeros(jkm)
  cp1_xyz = zeros(3); cp2_xyz = zeros(3); cp3_xyz = zeros(3)
  u1 = zeros(3); tvec = zeros(3)

  # loop over the CP planes in the direction di
  for jdi = 1:map.nctl[di]
    jkm[di] = jdi

    # set the (j,k,m) indices for the CPs at the corners of the plane
    cp1_jkm[di] = jdi
    cp1_jkm[it1] = 1
    cp1_jkm[it2] = 1

    cp2_jkm[di] = jdi
    cp2_jkm[it1] = 1
    cp2_jkm[it2] = map.nctl[it2]

    cp3_jkm[di] = jdi
    cp3_jkm[it1] = map.nctl[it1]
    cp3_jkm[it2] = 1

    # get the coordinates of the CPs at the corners of the plane, and some unit
    # vectors
    cp1_xyz = map.cp_xyz[:,cp1_jkm[1],cp1_jkm[2],cp1_jkm[3]]
    cp2_xyz = map.cp_xyz[:,cp2_jkm[1],cp2_jkm[2],cp2_jkm[3]]
    cp3_xyz = map.cp_xyz[:,cp3_jkm[1],cp3_jkm[2],cp3_jkm[3]]

    u1[:] = cp2_xyz - cp1_xyz
    fac1 = 1./norm(u1)
    u1[:] *= fac1

    tvec[it1] = -u1[it2]
    tvec[it2] = u1[it1]
    tvec[di] = 0.0
    alpha = dot(tvec, cp3_xyz - cp1_xyz)

    # u1 and tvec should be orthogonal
    @assert( dot(u1,tvec) < 1e-14 )

    # First equation, for coordinate it1
    # x3[it1] = -alpha*fac1*(x2[it2] - x1[it2]) + x1[it1]

    # set data for dependent cp3
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp3_idx[it1]
    LinG[ptr] = 1.0
    ptr += 1

    # set data for independent cp1
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp1_idx[it1]
    LinG[ptr] = -1.0
    ptr += 1
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp1_idx[it2]
    LinG[ptr] = -alpha*fac1
    ptr += 1

    # set data for independent cp2
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp2_idx[it2]
    LinG[ptr] = alpha*fac1
    ptr += 1

    Flow[fncidx] = 0.0
    Fupp[fncidx] = 0.0
    fncidx += 1

    # Second equation, for coordinate it2
    # x3[it2] = alpha*fac1*(x2[it1] - x1[it1]) + x1[it2]

    # set data for dependent cp3
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp3_idx[it2]
    LinG[ptr] = 1.0
    ptr += 1

    # set data for independent cp1
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp1_idx[it2]
    LinG[ptr] = -1.0
    ptr += 1
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp1_idx[it1]
    LinG[ptr] = alpha*fac1
    ptr += 1

    # set data for independent cp2
    iLfun[ptr] = fncidx
    jLvar[ptr] = cp2_idx[it1]
    LinG[ptr] = -alpha*fac1
    ptr += 1

    Flow[fncidx] = 0.0
    Fupp[fncidx] = 0.0
    fncidx += 1
  end

  return fncidx, ptr
end

@doc """
### numLinearStretchConstraints

Returns the number of stretch constraints (see setLinearStretchConstraints! for
further details).

**Arguments**

* `map` : PumiMapping object

"""->
function numLinearStretchConstraints(map::PumiMapping{Tffd}) where Tffd
  return map.nctl[1]*map.nctl[2]*map.nctl[3] - 1
end

@doc """
### countVarsLinearStretchConstraints!

Updates array `cntEq` with the number of variables involved in each
linear-stretch constraint.

**Arguments**

* `map` : PumiMapping object
* `cntEq` : array with number of variables involved in each equality constraint
* `fncidx` : index to the next equality constraint to add to cntEq

**Returns** `fncidx` + number of linear-plane constraints + 1

"""->
function countVarsLinearStretchConstraints!(map::PumiMapping{Tffd},
                                            cntEq::AbstractArray{Int,1},
                                            fncidx::Int) where Tffd
  num_cnstr = numLinearStretchConstraints(map)
  cntEq[fncidx:fncidx+num_cnstr-1] = 2
  return fncidx+num_cnstr
end

@doc """
### setLinearStretchConstraints!

This constraint is used to scale all the control points' `di` coordinate based
on a single master control point.  Useful when span is permitted to change and
you want each CP in a plane to stretch in proportion.

**Arguments**

* `map`  : PumiMapping object
* `di`   : the index the defines the constrained plane (e.g. `di`=2 -> j planes)
* `iLfun`, `jLvar` : the k^th element in the gradient matrix has row index
*                    `iLfun(k)` and column index `jLvar(k)`
* `LinG` : `LinG(k) = A(iLfun(k),jLvar(k))`
* `Flow`, `Fupp` : the upper and lower bounds on all constraints
* `fncidx` : function counter; next function to use in iLfun
* `ptr` : index pointer to next free space in iLfun and jLvar

**Returns** `fncidx` + number of linear-plane constraints + 1, and `ptr` + total
  number of non-zeros in the sparse contraints + 1

"""->
function setLinearStretchConstraints!(map::PumiMapping{Tffd}, di::Int, iLfun::AbstractArray{Int,1},
    jLvar::AbstractArray{Int,1}, LinG::AbstractArray{Tffd,1},
    Flow::AbstractArray{Tffd,1}, Fupp::AbstractArray{Tffd,1},
    fncidx::Int, ptr::Int) where Tffd

  # get master CP index and xyz
  master_jkm = ones(Int,3)
  master_jkm[di] = map.nctl[di]
  master_idx = map.cp_idx[di,master_jkm[1],master_jkm[2],master_jkm[3]]
  master_xyz = map.cp_xyz[di,master_jkm[1],master_jkm[2],master_jkm[3]]

  for k = 1:map.nctl[3]
    for j = 1:map.nctl[2]
      for i = 1:map.nctl[1]
        if i == master_jkm[1] && j == master_jkm[2] && k == master_jkm[3]
          # if this is the master node, skip it
          continue
        end

        # set data for dependent CP
        iLfun[ptr] = fncidx
        jLvar[ptr] = map.cp_idx[di,i,j,k]
        LinG[ptr] = 1.0
        ptr += 1

        # set data for master CP
        iLfun[ptr] = fncidx
        jLvar[ptr] = master_idx
        LinG[ptr] = map.cp_xyz[di,i,j,k]/master_xyz
        ptr += 1

        Flow[fncidx] = 0.0
        Fupp[fncidx] = 0.0
        fncidx += 1
      end
    end
  end
  return fncidx, ptr
end

@doc """
### numLinearRootConstraints

Returns the number of root constraints (see setLinearRootConstraints! for
further details).

**Arguments**

* `map` : PumiMapping object

"""->
function numLinearRootConstraints(map::PumiMapping{Tffd}) where Tffd
  return 5
end

@doc """
### countVarsLinearRootConstraints!

Updates array `cntEq` with the number of variables involved in each
linear-root constraint.

**Arguments**

* `map` : PumiMapping object
* `cntEq` : array with number of variables involved in each equality constraint
* `fncidx` : index to the next equality constraint to add to cntEq

**Returns** `fncidx` + number of linear-plane constraints + 1

"""->
function countVarsLinearRootConstraints!(map::PumiMapping{Tffd},
                                         cntEq::AbstractArray{Int,1},
                                         fncidx::Int) where Tffd
  for bdi = 1:3
    cntEq[fncidx] = 1
    fncidx += 1
  end
  cntEq[fncidx] = 2
  fncidx += 1
  cntEq[fncidx] = 2
  fncidx += 1
  return fncidx
end

@doc """
### setLinearRootConstraints!

This constraint sets some very specific linear constraints at the root of the
FFD, i.e. where the FFD touches the symmetric plane.  One constraint freezes one
CP completely, and another CP is keep aligned (in some direction) with the
frozen CP.  **Assumes the symmetry plane is at the low end of the CP indices**.

**Arguments**

* `map`  : PumiMapping object
* `di`   : defines the symmetry plane (e.g. `di`=2 -> y plane)
* `iLfun`, `jLvar` : the k^th element in the gradient matrix has row index
*                    `iLfun(k)` and column index `jLvar(k)`
* `LinG` : `LinG(k) = A(iLfun(k),jLvar(k))`
* `Flow`, `Fupp` : the upper and lower bounds on all constraints
* `fncidx` : function counter; next function to use in iLfun
* `ptr` : index pointer to next free space in iLfun and jLvar

**Returns** `fncidx` + number of linear-plane constraints + 1, and `ptr` + total
  number of non-zeros in the sparse contraints + 1

"""->
function setLinearRootConstraints!(map::PumiMapping{Tffd}, di::Int, iLfun::AbstractArray{Int,1},
    jLvar::AbstractArray{Int,1}, LinG::AbstractArray{Tffd,1},
    Flow::AbstractArray{Tffd,1}, Fupp::AbstractArray{Tffd,1},
    fncidx::Int, ptr::Int) where Tffd

  # frozen CP constraint
  for bdi = 1:3
    iLfun[ptr] = fncidx
    jLvar[ptr] = map.cp_idx[bdi,1,1,1]
    LinG[ptr] = 1.0
    ptr += 1

    Flow[fncidx] = map.cp_xyz[bdi,1,1,1]
    Fupp[fncidx] = Flow[fncidx]
    fncidx += 1
  end

  # Aligned CP: it is allowed to move in the direction it1 only
  it1 = mod(di,3)+1
  it2 = mod(di+1,3)+1
  jkm = ones(Int,3)
  jkm[it1] = map.nctl[it1]

  # constraint x2[di] = x1[di]
  iLfun[ptr] = fncidx
  jLvar[ptr] = map.cp_idx[di,jkm[1],jkm[2],jkm[3]]
  LinG[ptr] = 1.0
  ptr += 1
  iLfun[ptr] = fncidx
  jLvar[ptr] = map.cp_idx[di,1,1,1]
  LinG[ptr] = -1.0
  ptr += 1

  Flow[fncidx] = 0.0
  Fupp[fncidx] = 0.0
  fncidx += 1

  # constraint x2[it2] = x1[it2]
  iLfun[ptr] = fncidx
  jLvar[ptr] = map.cp_idx[it2,jkm[1],jkm[2],jkm[3]]
  LinG[ptr] = 1.0
  ptr += 1
  iLfun[ptr] = fncidx
  jLvar[ptr] = map.cp_idx[it2,1,1,1]
  LinG[ptr] = -1.0
  ptr += 1

  Flow[fncidx] = 0.0
  Fupp[fncidx] = 0.0
  fncidx += 1

  return fncidx, ptr
end
