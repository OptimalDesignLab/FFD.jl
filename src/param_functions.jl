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
  jkm = Array(Int, 3)
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


      end  # End for i = order+1:nctl
    end  # End for edg = 1:4

  end # End for di = 1:3

  return nothing
end  # end function calcChordLengthParam
