module FreeFormDeformation
export AbstractMappingType, Mapping, PumiMapping
export PumiBoundingBox, calcKnot, controlPoint, calcParametricMappingLinear
export calcParametricMappingNonlinear, evalVolume, findSpan

using ArrayViews
using PdePumiInterface
using ODLCommonTools
using SummationByParts


# Abstract Type definition
abstract AbstractMappingType # Abstract Mapping type for creating a different Mapping type
abstract AbstractBoundingBox

@doc """
### Mapping

A type for creating mapping objects needed for creating. The mapping is intended
for a uniform knot distribution along the 3 dimensions in the parametric space.
An object of this type can be defined by calling an inner constructor with the
following arguments in sequence

* Number of dimensions (2 or 3)
* Array of order of B-splines in all dimensions
* Number of control points along every direction
* Number of nodes (embedded geometry points) along every direction

**Method 2 constructor**
  function PumiMapping(ndim::Int, order::AbstractArray{Int,1},
                       nctl::AbstractArray{Int,1}, mesh_info::AbstractArray{Int,1};
                       full_geom=true, geom_faces=[0])

*  `ndim` : Number of dimensions
*  `order`: Aray of B-spline order in the 3 dimenstions
*  `nctl` : number of control points along the 3 dimenstions
*  `mesh_info` : information about the geometry being embedded. If the entire
                 geometry is being embedded then
```
  mesh_info = [num_nodes_per_element, num_elements]
```
                 else,
```
  mesh_info = [sbp.facenodes, num_boundary_faces, num_geometric_faces]
```
*  `full_geom` : Bool. `true` for embedding the entire geometry. `false` if you
                 want to embedded only certain faces of a geometry
*  `geom_faces` : Array of geometric face numbers being embedded in the FFD box

**Members**

*  `ndim`    : Number of dimensions (2D or 3D)
*  `nctl`    : Array of number of control points in each direction
*  `numnodes`: Array of number of nodes in each direction
*  `order`   : Array of order of B-splines in each direction
*  `xi`      : Array of parametric coordinates of the nodes of the geometry
               embedded within. The parametric axes within FFD are called
               (s,t,u)
*  `cp_xyz`  : Control point coordinates in the (x,y,z) space
*  `edge_knot` : 3D Array of edge knot vectors. dim1 = Knot vector, dim2 = edge
                 number (1:4), dim3 = direction (di)
*  `edge_param`: 3D array to store edge spacing parameters. dim1 = A or b,
                 dim2 = edge number (1:4), dim3 = direction (di)
*  `aj`   :
*  `dl` and `dr`   : Working array for computing spline basis and its derivatives
*  `work` :

"""->

type Mapping <: AbstractMappingType

  ndim::Int                     # Mapping object to indicate 2D or 3D
  nctl::AbstractArray{Int, 1}   # Number of control points in each of the 3 dimensions
  numnodes::AbstractArray{Int, 1} # Number of nodes in each direction
  order::AbstractArray{Int, 1}  # Order of B-spline in each direction

  xi::AbstractArray{AbstractFloat, 4}     # Coordinate values
  cp_xyz::AbstractArray{AbstractFloat, 4} # Cartesian coordinates of control points
  edge_knot::AbstractArray{Vector{AbstractFloat}, 1}  # edge knot vectors

  # Working arrays
  aj::AbstractArray{AbstractFloat, 3}
  dl::AbstractArray{AbstractFloat, 2}
  dr::AbstractArray{AbstractFloat, 2}
  work::AbstractArray{AbstractFloat, 4}

  evalVolume::Function

  function Mapping(dim, k, ncpts, nnodes)

    # Assertion statements to prevent errors
    @assert dim >= 2 "Only 2D and 3D valid"
    for i = 1:3
      @assert ncpts[i] > 0 "Number of control points specified is <= 0 in $i direction"
      @assert nnodes[i] > 0 "Number of nodes specified is <= 0 in $i direction"
      @assert k[i] > 0 "Order of B-spline specified is <= 0 in $i direction"
    end

    # Define max_wrk = number of work elements at each control point
    const max_work = 2*6  # 2*n_variables

    ndim = dim  # To indicate a 3D Mapping object is being created
    nctl = zeros(Int, 3)
    numnodes = zeros(Int, 3)
    order = zeros(Int, 3)

    nctl[:] = ncpts[:]    # Set number of control points
    numnodes[:] = nnodes[:] # Set map refinement level in each coordinate direction
    order[:] = k[:]       # Set the order of B-splines


    for i = 1:3
      @assert nctl[i] > 0
      @assert numnodes[i] > 0
      @assert order[i] > 0
    end

    # Allocate and initialize mapping arrays
    max_order = maximum(order)  # Highest order among 3 dimensions
    max_knot = max_order + maximum(nctl) # Maximum number of knots among 3 dimensions

    xi = zeros(numnodes[1], numnodes[2], numnodes[3], 3)
    cp_xyz = zeros(3,nctl[1], nctl[2], nctl[3])
    edge_knot = Array(Vector{AbstractFloat}, 3)
    for i = 1:3
      edge_knot[i] = zeros(AbstractFloat, nctl[i]+order[i])
    end

    aj = zeros(3, max_order, 3)
    dl = zeros(max_order-1, 3)
    dr = zeros(max_order-1, 3)
    work = zeros(nctl[1], nctl[2], nctl[3], max_work)

    new(ndim, nctl, numnodes, order, xi, cp_xyz, edge_knot, aj, dl, dr, work)

  end  # End constructor

end  # End Mapping

type PumiMapping{Tffd} <: AbstractMappingType

  ndim::Int                     # Mapping object to indicate 2D or 3D
  full_geom::Bool               # Embed entire geometry or only certain faces
  nctl::AbstractArray{Int, 1}   # Number of control points in each of the 3 dimensions
  order::AbstractArray{Int, 1}  # Order of B-spline in each direction

  xi::AbstractArray{Tffd}       # Paramaetric coordinates of input geometry
  cp_xyz::AbstractArray{Tffd, 4} # Cartesian coordinates of control points
  edge_knot::AbstractArray{Vector{Tffd}, 1}  # edge knot vectors
  geom_faces::AbstractArray{Int,1}

  # Working arrays
  aj::AbstractArray{Tffd, 3}
  dl::AbstractArray{Tffd, 2}
  dr::AbstractArray{Tffd, 2}
  work::AbstractArray{Tffd, 4}

  evalVolume::Function

  function PumiMapping(ndim::Int, order::AbstractArray{Int,1},
                       nctl::AbstractArray{Int,1}, mesh_info::AbstractArray{Int,1};
                       full_geom=true, geom_faces::AbstractArray{Int,1}=[0])

    map = new()
    # Check if the input arguments are valid
    @assert ndim >= 2 "Only 2D and 3D valid"
    for i = 1:3
      @assert order[i] > 0 "Order cannot be 0 in the $i direction"
      @assert nctl[i] > 0 "number of control points cannot be 0 in $i direction"
    end
    @assert mesh_info[1] > 0 "Number of nodes per element can only be > 0"
    @assert mesh_info[2] > 0 "Grid to be deofrmed must have > 0 elements"

    # Define max_wrk = number of work elements at each control point
    const max_work = 2*6  # 2*n_variables

    map.ndim = ndim
    map.full_geom = full_geom
    map.order = order
    map.nctl = nctl
    map.cp_xyz = zeros(3, nctl[1], nctl[2], nctl[3])

    map.edge_knot = Array(Vector{Tffd}, 3)
    for i = 1:3
      map.edge_knot[i] = zeros(Tffd, nctl[i]+order[i])
    end

    if full_geom == true
      map.xi = zeros(Tffd, 3, mesh_info[1], mesh_info[2])
    else
      map.xi = zeros(Tffd, 3, mesh_info[1], mesh_info[2], mesh_info[3])
      map.geom_faces = geom_faces
    end

    # Allocate and initialize mapping arrays
    max_order = maximum(order)  # Highest order among 3 dimensions
    max_knot = max_order + maximum(nctl) # Maximum number of knots among 3 dimensions

    map.aj = zeros(3, max_order, 3)
    map.dl = zeros(max_order-1, 3)
    map.dr = zeros(max_order-1, 3)
    map.work = zeros(nctl[1], nctl[2], nctl[3], max_work)

    return map
  end
end  # End type PumiMapping

include("knot.jl")
include("bounding_box.jl")
include("mapping_functions.jl")
include("control_point.jl")
include("span.jl")
include("b-splines.jl")
include("evaluations.jl")

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
    calcKnot(map)

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
          map.aj[:,jc3,3] = map.cp_xyz[p,q,imk[3]+jc3,:]
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

function contractWithdGdB(map, dJdGrid)

  @assert map.numnodes[1] == size(dJdGrid, 1)
  @assert map.numnodes[2] == size(dJdGrid, 2)
  @assert map.numnodes[3] == size(dJdGrid, 3)

  # Allocate and initialize arrays
  basis = zeros(AbstractFloat, maximum(map.order), 3)
  fill!(map.work, 0.0)

  # Loop over the nodes of the mapping
  for m = 1:map.numnodes[3]
    for k = 1:map.numnodes[2]
      for j = 1:map.numnodes[1]

        # Store the parameter values and dJdGrid
        xi[:] = map.xi[j,k,m,:]
        dJdG = view(dJdGrid, j, k, m, 1:3)

        # Evaluate the knot vectors and basis values
        # TODO: Compute the knot vector for subsequent operations
        for di = 1:3
          span[di] = findSpan(xi[di], map, di)
          basisFunctions(map, basis[:,di], di, xi[di], span[di])
        end  # End for di = 1:3

        for r = 1:map.order[3]
          rs = r + span[3] - map.order[3]
          for q = 1:map.order[2]
            qs = q + span[2] - map.order[2]
            for p = 1:map.order[1]

              ps = p + span[1] - map.order[1]
              coeff = basis[p,1]*basis[q,2]*basis[r,3]
              map.work[ps, qs, rs, 1:3] += coedd*dJdG[:]

            end # End for p = 1:map.order[1]
          end   # End for q = 1:map.order[2]
        end     # End for r = 1:map.order[3]

      end  # End for j = 1:map.numnodes[1]
    end    # End for k = 1:map.numnodes[2]
  end      # End for m = 1:map.numnodes[3]

  return nothing
end  # End function contractWithdGdB(map, dJdGrid)

end # End module FreeFormDeformation
