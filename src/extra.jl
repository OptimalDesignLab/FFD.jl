# Superflous code
# Contains code ported from Prof. Hicksn's Mapping_Mod.f90

@doc """
### calcJacobian

!!!!!   NOT NEEDED   !!!!
Calculates the metric Jacobian value based on the B-spline explicit mapping.

**Arguments**

*  `map`  : Object of Mapping type
*  `Jac`  : Array of mesh Jacobians in j,k,m format

"""->

function calcJacobian(map, Jac)

  # Check if map.jkmmax = size(Jac)
  @assert map.jkmmax[1] == size(Jac,1)
  @assert map.jkmmax[2] == size(Jac,2)
  @assert map.jkmmax[3] == size(Jac,3)

  # Allocate and initialize arrays
  scal = zeros(AbstractFloat, 3)
  xi = zeros(AbstractFloat, 3)
  dx = zeros(AbstractFloat, 3, 3)

  scal[:] = 1 ./ (map.jkmmax[:] - 1)

  for m = 1:map.jkmmax[3]
    for k = 1:map.jkmmax[2]
      for j = 1:map.jkmmax[1]
        xi[:] = map.jkmmax[j,k,m,:]
        calcdXdxi(map, xi, jderiv, dx[:,1])
        calcdXdxi(map, xi, jderiv, dx[:,2])
        calcdXdxi(map, xi, jderiv, dx[:,3])
        dx[:,1] *= scal[1]
        dx[:,2] *= scal[2]
        dx[:,3] *= scal[3]
        Jac[j,k,m] = dx[1,1]*dx[2,2]*dx[3,3] + dx[1,2]*dx[2,3]*dx[3,1] +
                     dx[1,3]*dx[2,1]*dx[3,2] - dx[1,1]*dx[2,3]*dx[3,2] -
                     dx[1,2]*dx[2,1]*dx[3,3] - dx[1,3]*dx[2,2]*dx[3,1]
        Jac[j,k,m] = 1/Jac[j,k,m]
      end  # End for j = 1:map.jkmmax[1]
    end  # End for k = 1:map.jkmmax[2]
  end  # End for m = 1:map.jkmmax[3]

  return nothing
end

@doc """
### resizeMapping

!!!! NOT NEEDED FOR FFD  !!!!

Uses edge spacing parameters to refine/coarsen in each parameter direction.
This might be useful for refinement studies where we want to insert or delete
points in a continous way, or for grid sequencing/multigrid applications.

**Arguments**

*  `map`  : Object of Mapping type
*  `jkmmax` : dersired number of (xi, eta, zeta) parameters

"""->

function resizeMapping(map, jkmmax)

  for i = 1:3
    @assert jkmmax[i] > 0
  end

  # Check if resizing is necessary
  if map.jkmmax[1] != jkmmax[1] || map.jkmmax[2] != jkmmax[2] || map.jkmmax[3]
     != jkmmax[3]
    map.xi = zeros(AbstractFloat, jkmmax[1], jkmmax[2], jkmmax[3], 3)
  end

  uvw = zeros(AbstractFloat, 3)

  # Loop over all nodes, calculate each new (xi, eta, zeta)
  for m = 1:map.jkmmax[3]
    uvw[3] = (m-1)/(map.jkmmax[3] - 1)
    for k = 1:map.jkmmax[2]
      uvw[2] = (k-1)/(map.jkmmax[2]-1)
      for j = 1:map.jkmmax[1]

        uvw[1] = (j-1)/(map.jkmmax[1] - 1)
        calcXi(map, uvw, map.xi[j,k,m,:])
        # need to set lower/upper value explicitly to avoid round-off errors
        if j == 1
          map.xi[j,k,m,1] = 0.0
        end
        if k == 1
          map.xi[j,k,m,2] = 0.0
        end
        if m == 1
          map.xi[j,k,m,3] = 0.0
        end
        if j == map.jkmmax[1]
          map.xi[j,k,m,1] = 1.0
        end
        if k == map.jkmmax[2]
          map.xi[j,k,m,2] = 1.0
        end
        if m == map.jkmmax[3]
          map.xi[j,k,m,3] = 1.0
        end

      end  # End for j = 1:map.jkmmax[1]
    end  # End for k = 1:map.jkmmax[2]
  end  # End for m = 1:map.jkmmax[3]

  return nothing
end  # End function resizeMapping(map, jkmmax)

@doc """
### calcXi

!!!! NOT NEEDED FOR FFD !!!

Calculates the intermediate parameters, (xi,eta,zeta), using the edge spacing
control functions defined by map%edgpar(1:2,:,:)

**Arguments**

*  `map`  : Object of Mapping type
*  `uvw`  : j,k,m indices scaled to fit in [0,1],  i.e. uniform spacing
*  `xi`   : Parametric coordinates (length = 3)

"""->

function calcXi(map, uvw, xi)

  # Allocate and initialize arrays
  # xi = zeros(AbstractFloat, 3)
  xiedg = zeros(AbstractFloat, 4, 3)
  F = zeros(AbstractFloat, 3)
  dF = zeros(AbstractFloat, 3, 3)
  dxi = zeros(AbstractFloat, 3)
  psi = zeros(AbstractFloat, 3)
  dpsi = zeros(AbstractFloat, 3)


  for di = 1:3
    it1 = mod(di,3) + 1
    it2 = mod(di+1,3) + 1
    for edg 1:4
      A = map.edge_param[1, edg, di]
      b = map.edge_param[2, edg, di]
      fac = uvw[di] - 0.5
      if b == 0.0
        # Use taylor series to determine u
        fac = 1.0 + 2.0*fac
      else
        fac = 1.0 + tanh(b*fac)/tanh(0.5*b)
      end  # End if b == 0.0
      xiedg[edg, di] = fac/( 2*A + (1 - A)*fac )
    end  # End for edg 1:4
  end  # End for di = 1:3

  # calculate the intermediate params by find the intersection of 3 trilinear
  # surface (via Newton's method)

  # Get initial estimate
  for di = 1:3
    xi[di] = 0.0
    for edg = 1:4
      xi[di] = 0.25*xiedg[edg,di]
    end
  end

  # loop until sufficiently converged

  for n = 1:100
    # Get the RHS and LHS for Newton's method
    F[:] = -xi[:]
    fill!(dF, 0.0)
    for di = 1:3
      it1 = mod(di,3) + 1
      it2 = mod(di+1,3) + 1
      dF[di,di] = 1.0
      for edg = 1:4
        if mod(edg,2) == 1
          psi[it1] = 1.0 - xi[it1]
          dpsi[it1] = 1.0
        else
          psi[it1] = xi[it1]
          dpsi[it1] = 1.0
        end  # End if mod(edg,2) == 1

        if mod((edg+1)/2,2) == 1
          psi[it2] = 1.0 - xi[it2]
          dpsi[it2] = -1.0
        else
          psi[it2] = xi[it2]
          dpsi[it2] = 1.0
        end  # End if mod((edg+1)/2,2) == 1

        F[di] += psi[it1]*psi[it2]*xiedg[edg,di]
        dF[di,it1] -= -dpsi[it1]*psi[it2]*xiedg[edg,di]
        dF[di,it2] -= -psi[it1]*dpsi[it2]*xiedg[edg,di]
      end  # End for edg = 1:4
    end # End for di = 1:3

    # Solve for the update
    # TODO: Solve for the update

    # update parameter and exit loop
    xi[:] += F[:]
    if norm(F,2) < 1e-13
      return nothing
    end

  end  # End for n = 1:100

  return nothing
end

# parameter functions
@doc """
### calcChordLengthParam

!!! NOT NEEDED FOR FFD !!!

Calculates the mapping parameter values based on the chord lengths in physical
xyz space. It also calculates the edge knot vectors.

**Inputs**

*  `map` : Object of Mapping type
*  `xyz` : Mesh in physical space on which the chord length parameterization is
           based on

**Outputs**

*  None

""" ->

function calcChordLengthParam(map, xyz)

  # Intermediate variables
  jkm = zeros(Int, 3)
  vec = zeros(3)

  # Loop over each direction, Find the parameter values for a curve of constant
  # it1 & it2
  for di = 1:3
    it1 = mod(di,3) + 1
    it2 = mod(di+1,3) + 1

    for jit1 = 1:jkmmax[it1]
      jkm[it1] = jit1
      for jit2 = 1:jkmmax[it2]
        jkm[it2] = jit2
        jkm[di] = 1
        map.xi[jkm[1], jkm[2], jkm[3], di] = 0.0
        for jdi = 1:(map.jkmmax[di]-1)
          jkm[di] = jdi
          j = jkm[1]
          k = jkm[2]
          m = jkm[3]
          jkm[di] = jdi + 1
          jp = jkm[1]
          kp = jkm[2]
          mp = jkm[3]

          # Get the arc-length element
          vec[:] = xyz[jp, kp, mp, :] - xyz[j, k, m, :]
          ds = norm(vec,2)
          map.xi[jp,kp,mp,di] = map.xi[j,k,m,di] + ds
        end  # End for jdi = 1:(map.jkmmax[di]-1)

        # Normalize the arc-length parameters
        jkm[di] = map.jkmmax[di]
        j = jkm[1]
        k = jkm[2]
        m = jkm[3]
        fac = 1 / map.xi[j,k,m,di]
        map.xi[j,k,m,di] = 1.0
        for jdi = 2:(map.jkmmax[di]-1)
          jkm[di] = jdi
          j = jkm[1]
          k = jkm[2]
          m = jkm[3]
          map.xi[j,k,m,di] = fac*map.xi[j,k,m,di]
        end  # End for jdi = 2:(map.jkmmax[di]-1)
      end  # End for jit2 = 1:jkmmax[it2]
    end  # End for jit1 = 1:jkmmax[it1]
  end  # end for di = 1:3

  # At this point, we have arc length parameters, so arc-length-based edge knots
  # can be found

  for di = 1:3
    it1 = mod(di,3) + 1
    it2 = mod(di+1,3) + 1
    order = map.order[di]
    nctl = map.nctl[di]

    for edg = 1:4
      if mod(edg,2) == 1
        jkm[it1] = 1
      else
        jkm[it1] = map.jkmmax[it1]
      end  # End if mod(edg,2) == 1

      if mod((edg+1)/2,2) == 1
        jkm[it2] = 1
      else
        jkm[it2] = map.jkmmax[it2]
      end  # End if mod((edg+1)/2,2) == 1

      # Set the edge knot vector's end knots
      map.edge_knot[1:order, edg, di] = 0.0
      map.edge_knot[nctl+1:nctl+order, edg, di] = 1.0

      # Find point intervals per knot intervals
      ppk = (map.jkmmax[di]-1) / (nctl+1-order)

      for i = order+1:nctl
        knotloc = ppk*(i-order)  # knot location
        for jdi = 1:(map.jkmmax[di]-1)
          jkm[di] = jdi

          if knotloc >= jdi-1 && knotloc < jdi
            jkm[di] = jdi
            j = jkm[1]
            k = jkm[2]
            m = jkm[3]
            jkm[di] = jdi + 1
            jp = jkm[1]
            kp = jkm[2]
            mp = jkm[3]
            map.edge_knot[i,edg,di] = (map.xi[jp,kp,mp,di] - map.xi[j,k,m,di])*
                                      (knotloc - jdi + 1) + map.xi[j,k,m,di]
          end  # End if knotloc >= jdi-1 && knotloc < jdi
        end # End for jdi = = 1:(map.jkmmax[di]-1)
      end  # End for i = order+1:nctl
    end  # End for edg = 1:4
  end # End for di = 1:3

  return nothing
end  # end function calcChordLengthParam

@doc """
### calcEdgeSpacingParam


!!!  NOT NEEDED  !!!
Uses the xi, eta, zeta values to calculate the parameters for hyperbolic tangent
spacing function. These parameters are

A = sqrt(sp2/sp1)
b, where, sinh(b) - b/((N-1) \* sp1 \* sp2) = 0

and sp1 and sp2 are the (xi, eta, zeta) spacing at the ends of the edges.
A is stored in `map.edge_param[1,:,:]`, and b in `map.edge_param[2, :, :]`

**Arguments**

*  `map` : Object of Mapping type

""" ->

function calcEdgeSpacingParam(map)

  # Initialize working array
  jkm = zeros(Int,3)

  for di = 1:3
    it1 = mod(di,3) + 1
    it2 = mod(di+1,3) + 1

    for edg = 1:4
      if mod(edg,2) == 1
        jkm[it1] = 1
      else
        jkm[it1] = map.jkmmax[it1]
      end

      if mod((edg+1)/2,2) == 1
        jkm[it2] = 1
      else
        jkm[it2] = map.jkmmax[it2]
      end

      # Find the spacing at the ends of the edge edg
      jkm[di] = 1
      j = jkm[1]
      k = jkm[2]
      m = jkm[3]
      jkm[di] = 2
      jp = jkm[1]
      kp = jkm[2]
      mp = jkm[3]
      dx2 = map.xi[jp,kp,mp,di] - map.xi[j,k,m,di]

      ### #if 0 satement is supposed to come here
      ###
      ###

      # Set parameter A
      map.edge_param[1, edg,di] = sqrt(dx2/dx1)

      # Set parameter b
      Nm1 = map.jkmmax[di] - 1
      x2 = (map.jkmmax[di] - 1)*sqrt(dx1*dx2)
      if x2 > 1
        # 1/(N-1) is smaller than sqrt(dx1*dx2), so b will be zero
        map.edge_param[2, edg, di] = 0.0
      else
        # need to find b using Newton's method
        x1 = small
        x2 = sqrt(6*(1/x2 - 1))
        map.edge_param[2,edg,di] =  bfuncNewtonSolve(x2, Nm1, dx1, dx2)
      end  # End if x2 > 1
    end    # End for edg = 1:4
  end      # End for di = 1:3

  return nothing
end  # End function calcEdgeSpacingParam(map)

function bfunc(x, Nm1, dx1, dx2)

  tmp = 1 / (Nm1*sqrt(dx1*dx2))
  f = sinh(x) - x*tmp
  df = cosh(x) - tmp

  return f, df
end  # end function bfunc

function bfuncNewtonSolve(guess, Nm1, dx1, dx2)

  b = guess
  b_next = zero(AbstractFloat)

  for i = 1:50
    f, df = bfunc(b, Nm1, dx1, dx2)
    b_next = b - f/df
    if norm(b_next-b, 2) < 1e-12
      break
    else
      b = b_next
    end
  end # End for i = 1:50

  return b_next
end


# Grid.jl

####  -------------------

#  NOT NEEDED !!!!!!!!!!

#### ---------------------

@doc """
### calcGrid

Calculates the mesh coordinates based on the element node values in the map,
hermite basis functions and the refinement level. The basis function values of
map should be known

**Arguments**

*  `map` : Object of Mapping type
*  `xyz` : Array of mesh node coordinates in j,k,m,1:3 format

"""->

function calcGrid(map, xyz)

  @assert size(map.jkmmax) == size(xyz)

  for m = 1:map.jkmmax[3]
    for k = 1:map.jkmmax[2]
      for j = 1:map.jkmmax[1]
        xi = view(map.ki, j, k, m, :)
        dXdxi(map, xi, jderiv, xyz[j,k,m,:])
      end  # End for j = 1:map.jkmmax[1]
    end  # End for k = 1:map.jkmmax[2]
  end  # End for m = 1:map.jkmmax[3]

  return nothing
end


@doc """
### fitGrid

Finds the control point values such that, when evaluated at in computational
space, the mapping produces the mesh xyz as close as possible. Parameter
correction was added by L. Olague and J. Hicken and ported to the language
`Julia`

**Arguments**

*  `map` : Object of Mapping type
*  `xyz` : Mesh to be approximated. On exit, htis is over written with the new
           mesh produced by the map
*  `L2_error` : The L2 error squared in fitting the grid is added to this
                variable at the end of routine
*  `max_error`: The values in this array are updated if the infinity norm in the
                fit of this grid is larger

"""


function fitGrid(map, xyz, L2_error, max_error)

  # Allocate and initialize arrays


  return nothing
end
