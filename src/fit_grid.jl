# fitGrid()
@doc """
### fitGrid

Finds the control point values such that, when evaluated at in computational
space, the mapping produces the mesh xyz as close as possible. Does not have any
parameter correction currently

**Arguments**

*  `map` : Object of Mapping type
*  `xyz` : The mesh to be approximated. Upon exit, it is overwritten with the
           new mesh produced by the map.

"""->

function fitGrid(map, xyz)

  # Get initial chord length parameterization and edge knot vectors
  calcChordLengthParam(map, xyz)

  # Maximum number of control points being solved for
  N = (map.nctl[1] - 2)*(map.nctl[2] - 2)*(map.nctl[3] - 2)
  # Maximum number of nodes in the RHS vector
  Nnodes = map.jkmmax[1]*map.jkmmax[2]*map.jkmmax[3]
  # Number of control points that influence a node
  maxcols = map.order[1]*map.order[2]*map.order[3]

  # Allocate and initialize arrays
  index = zeros(Int, map.nctl[1], map.nctl[2], map.nctl[3])
  ia = zeros(Int, Nnodes+1)
  ja = zeros(Int, Nnodes*maxcols)
  aa = zeros(AbstractFloat, Nnodes*maxcols)
  iat = zeros(Int, N+1)
  jat = zeros(Int, Nnodes*maxcols)
  aat = zeros(AbstractFloat, Nnodes*maxcols)
  rhs = zeros(AbstractFloat, Nnodes, 3)
  S = zeros(AbstractFloat, N, N)
  b = zeros(AbstractFloat, N, 3)
  Ncoeff = zeros(AbstractFloat, maximum(map.jkmmax), maximum(map.nctl))
  basis = zeros(AbstractFloat, maximum(map.order), 3)
  jkm = zeros(Int, 3)
  jkmcp = zeros(Int, 3)
  span = zeros(Int, 3)
  xi = zeros(AbstractFloat, 3)

  # Set parameter spacing
  dxi = zeros(Int, 3)
  dxi = 1./(map.jkmmax - 1)

  # initialize
  fill!(map.cp_xyz, 0.0)

  # The edges should be fit first to ensure consistency between blocks that
  # share edges. Next faces should be fitted followed by the volumes

  # Loop over each edge
  for di = 1:3
    it1 = mod(di,3) + 1
    it2 = mod(di+1,3) + 1
    for edg = 0:3
      if mod(edg,2) == 0
        jkm[it1] = 1
        xi[it1] = 0.0
      else
        jkm[it1] = map.jkmmax[it1]
        xi[it1] = 1.0
      end  # End if mod(edg,2) == 0

      if mod(edg/2,2) == 0
        jkm[it2] = 1
        xi[it2] = 0.0
      else
        jkm[it2] = map.jkmmax[it2]
        xi[it2] = 1.0
      end  # End if mod(edg/2,2) == 0

      # Form the matrix, N, of B-spline basis functions
      fill!(Ncoeff, 0.0)

      for jdi = 1:map.jkmmax(di)
        jkm[di] = jdi
        j = jkm[1]
        k = jkm[2]
        m = jkm[3]

        xi[di] = map.xi[j, k, m, di] # get the parameter value

        # Get the Control Point coefficients
        calcKnot(map, xi)  # TODO: Create function called calcKnot
        span[di] = findSpan(xi[di], map, di)
        basisFunctions(map, view(basis,:,di), di, xi[di], span[di])
        for i = 1:map.order[di]
          Ncoeff[jdi, i+span[di]-map.order[di]] = basis[i,di]
        end

        rhs[jdi,:] = xyz[j,k,m,:]  # Set the RHS value
      end  # End for jdi = 1:map.jkmmax(di)

      # Solve for the control points that minimizes the least squares problem

    end  # End for edg = 0:3

  end  # End for di = 1:3

  return nothing
end
