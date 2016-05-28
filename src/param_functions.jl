# parameter functions
@doc """
### calcChordLengthParam

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
        map.edge_param[2,edg,di] = ### TODO: figure out newton's method
      end  # End if x2 > 1
    end    # End for edg = 1:4
  end      # End for di = 1:3

  return nothing
end

function bfunc(x, f, df, Nm1, dx1, dx2)

  tmp = 1 / (Nm1*sqrt(dx1*dx2))
  f = sinh(x) - x*tmp
  df = cosh(x) - tmp

  return f, df
end  # end function bfunc
