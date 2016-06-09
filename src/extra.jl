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
