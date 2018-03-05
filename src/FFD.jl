module FFD

using MPI
# initialize MPI, and arrange for its finalization if so
function finalizeMPI()
#  println("running atexit hook for MPI")
  if MPI.Initialized()
#    println("finalizing MPI")
    MPI.Finalize()
  end
end

if !MPI.Initialized()
  println("initialiing MPI")
  MPI.Init()
  atexit(finalizeMPI)
end


export AbstractMappingType, Mapping, PumiMapping
export PumiBoundingBox, calcKnot, controlPoint, calcParametricMappingLinear
export calcParametricMappingNonlinear, evalVolume, evalSurface
export writeControlPointsVTS, evaldXdControlPointProduct
export initializeFFD, commitToPumi

export numLinearPlaneConstraints, countVarsLinearPlaneConstraints!, setLinearPlaneConstraints!
export numLinearCornerConstraints, countVarsLinearCornerConstraints!, setLinearCornerConstraints!
export numLinearStretchConstraints, countVarsLinearStretchConstraints!, setLinearStretchConstraints!
export numLinearRootConstraints, countVarsLinearRootConstraints!, setLinearRootConstraints!

push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))

using ArrayViews
using PdePumiInterface
using ODLCommonTools
using SummationByParts
using WriteVTK


# Abstract Type definition
abstract AbstractMappingType{Tffd} <: Any # Abstract Mapping type for creating a different Mapping type
abstract AbstractBoundingBox{Tffd} <: Any


@doc """
### Mapping

A type for creating mapping objects. The mapping uses a uniform knot
distribution along the 3 dimensions in the parametric space.  An object of this
type can be defined by calling an inner constructor with the following arguments
in sequence

* Number of dimensions (2 or 3)
* Array of order of B-splines in all dimensions
* Number of control points along every direction
* Number of nodes (embedded geometry points) along every direction

**Method 2 constructor**
  function PumiMapping(ndim::Int, order::AbstractArray{Int,1},
                       nctl::AbstractArray{Int,1}, mesh_info::AbstractArray{Int,1};
                       full_geom=true, bc_nums=[0])

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
*  `bc_nums` : Array of boundary condition numbers that identify the surface 
               being embedded in the FFD box

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

#TODO: get rid of abstract types
@doc """
### PumiMapping

"""->

type PumiMapping{Tffd} <: AbstractMappingType{Tffd}

  ndim::Int                     # Mapping object to indicate 2D or 3D
  full_geom::Bool               # Embed entire geometry or only certain faces
  nctl::AbstractArray{Int, 1}   # Number of control points in each of the 3 dimensions
  order::AbstractArray{Int, 1}  # Order of B-spline in each direction

  xi::AbstractArray       # Paramaetric coordinates of input geometry
  cp_xyz::AbstractArray{Tffd, 4} # Cartesian coordinates of control points
  edge_knot::AbstractArray{Vector{Tffd}, 1}  # edge knot vectors
  bc_nums::AbstractArray{Int,1}
  cp_idx::AbstractArray{Int,4}  # index assigned to each CP coordinate

  # Working arrays
  aj::AbstractArray{Tffd, 3}
  dl::AbstractArray{Tffd, 2}
  dr::AbstractArray{Tffd, 2}
  work::AbstractArray{Tffd, 4}

  evalVolume::Function

  function PumiMapping(ndim::Int, order::AbstractArray{Int,1},
                       nctl::AbstractArray{Int,1}, mesh::AbstractMesh;
                       full_geom=true, bc_nums::AbstractArray{Int,1}=[0])

    map = new()
    # Check if the input arguments are valid
    @assert ndim >= 2 "Only 2D and 3D valid"
    @assert mesh.coord_order == 1
    for i = 1:3
      @assert order[i] > 0 "Order cannot be 0 in the $i direction"
      @assert nctl[i] > 0 "number of control points cannot be 0 in $i direction"
    end

    # Define max_wrk = number of work elements at each control point
    const max_work = 2*6  # 2*n_variables

    map.ndim = mesh.dim
    map.full_geom = full_geom
    map.order = order
    map.nctl = nctl
    map.cp_xyz = zeros(3, nctl[1], nctl[2], nctl[3])

    map.edge_knot = Array(Vector{Tffd}, 3)
    for i = 1:3
      map.edge_knot[i] = zeros(Tffd, nctl[i]+order[i])
    end

    if full_geom == true
      map.xi = zeros(Tffd, 3, mesh.dim+1, mesh.numEl) # Only valid for simplex elements
    else
      map.bc_nums = bc_nums
      map.xi = Array(Array{Tffd,3}, length(bc_nums))
      defineMapXi(mesh, bc_nums, map.xi)
    end

    # Allocate and initialize mapping arrays
    max_order = maximum(order)  # Highest order among 3 dimensions
    max_knot = max_order + maximum(nctl) # Maximum number of knots among 3 dimensions

    # use simple logical indexing for CP coordinates
    map.cp_idx = zeros(3, nctl[1], nctl[2], nctl[3])
    ptr = 1
    for k = 1:nctl[3]
      for j = 1:nctl[2]
        for i = 1:nctl[1]
          for di = 1:3
            map.cp_idx[di,i,j,k] = ptr
            ptr += 1
          end
        end
      end
    end

    map.aj = zeros(3, max_order, 3)
    map.dl = zeros(max_order-1, 3)
    map.dr = zeros(max_order-1, 3)
    map.work = zeros(max_work, nctl[1], nctl[2], nctl[3])

    return map
  end

  function PumiMapping()
    ndim = 0
    full_geom = false
    nctl = Int[]
    order = Int[]
    cp_xyz = zeros(Tffd, 0, 0, 0, 0)
    edge_knot = Array(Vector{TFFD}, 0)
    bc_nums = Int[]
    cp_idx = zeros(Tffd, 0, 0, 0, 0)
    aj = zeros(Tffd, 0, 0, 0)
    dl = zeros(Tffd, 0, 0)
    dr = zeros(Tffd, 0, 0)
    work = zeros(Tffd, 0, 0, 0, 0)
    evalVolume = () -> nothing
    

    return new(ndim, full_geom, nctl, order, cp_xyz, edge_know, bc_nums,
               cp_idx, aj, dl, dr, work, evalVolume)
end



end  # End type PumiMapping



include("knot.jl")
include("bounding_box.jl")
include("mapping_functions.jl")
include("control_point.jl")
include("span.jl")
include("b-splines.jl")
include("evaluations.jl")
include("constraints.jl")
include("pumi_specific_functions.jl")

@doc """
Routine to be called externally for initializing FreeFormDeformation

**Input**

* `mesh` : Pumi mesh
* `sbp`  : Summation-By-Parts operator
* `order`: Array of B-splone order along the 3 parametric directions (ξ, η, ζ)
* `nControlPts` : Array of nu,ber of control points along 3 parametric
                  directions (ξ, η, ζ)
* `offset` : FFD Bounding box offset form the embedded geometry along the 3
             physical directions (x, y, z)
* `full_geom` : Whether the full geometry or a portion of geometry is being
                embedded. (True or false)
* `bc_nums`: (Optional Argument) If partial geometry embedded, give the
                face/edge numbers of the embedded gemetry

"""->

function initializeFFD{Tmsh}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                       order::AbstractArray{Int,1},
                       nControlPts::AbstractArray{Int,1},
                       offset::AbstractArray{Float64,1}, full_geom::Bool,
                       bc_nums::AbstractArray{Int,1}=[0])

  # Create Mapping object
  ndim = mesh.dim
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=full_geom,
                              bc_nums=bc_nums)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  ffd_box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, ffd_box)

  # Populate map.xi
  if full_geom == true
    calcParametricMappingNonlinear(map, ffd_box, mesh)
  else
    calcParametricMappingNonlinear(map, ffd_box, mesh, bc_nums)
  end

  return map, ffd_box
end

@doc """
### FreeFormDeformation.defineMapXi

"""->

function defineMapXi(mesh::AbstractMesh, bc_nums::AbstractArray{Int,1},
                     xi::AbstractArray)

  for (idx, itr) in enumerate(bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = view(mesh.bndryfaces, idx_range)
    nfaces = length(bndry_facenums)
    xi[idx] = zeros(3,mesh.coord_numNodesPerFace,nfaces)
  end

  return nothing
end

@doc """
###FreeFormDeformation.defineVertices

This function defines the shape of the vertices array that is used to update
the Pumi mesh after FFD when a geometric face is paramtereized.

**Arguments**

* `mesh` : Pumi DG mesh
* `bc_nums` : Array of boundary condition numbers over which the number
                 of element faces needs to be computed
* `vertices` : Array of arrays holding the coodinates of the updated vertices
               shape = vertices[n_bc_nums][mesh.dim, mesh.coord_numNodesPerFace, n_elem_faces]

"""->

function defineVertices{Tmsh}(mesh::AbstractDGMesh{Tmsh}, bc_nums::AbstractArray{Int,1},
                        vertices::AbstractArray)

  for (idx, itr) in enumerate(bc_nums)
    start_index = mesh.bndry_offsets[itr]
    end_index = mesh.bndry_offsets[itr+1]
    idx_range = start_index:(end_index-1)
    bndry_facenums = sview(mesh.bndryfaces, idx_range)
    nfaces = length(bndry_facenums)
    arr_dim = mesh.dim
    vertices[idx] = zeros(Tmsh, arr_dim, mesh.coord_numNodesPerFace, nfaces)
    for i = 1:nfaces
      bndry_i = bndry_facenums[i]
      for j=1:mesh.coord_numNodesPerFace
        v_j = mesh.topo.face_verts[j, bndry_i.face]
        vertices[idx][1:mesh.dim,j,i] = mesh.vert_coords[1:mesh.dim,v_j,bndry_i.element]
      end
    end # End for i = 1:nfaces
  end

  return nothing
end

#=
@doc """
###FreeFormDeformation.commitToPumi

Uses the modified vertex coordinates to update the Pumi mesh

"""->

function commitToPumi{Tffd, Tmsh}(map::PumiMapping{Tffd},
                      mesh::AbstractDGMesh{Tmsh}, sbp, vertices, opts)

  if map.full_geom == true
    println("Tmsh = $Tmsh")
    for i = 1:mesh.numEl
      update_coords(mesh, i, real(vertices[:,:,i]))
    end
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
      idx_range = start_index:(end_index-1)
      bndry_facenums = view(mesh.bndryfaces, idx_range) # faces on geometric edge i
      nfaces = length(bndry_facenums)
      for i = 1:nfaces
        bndry_i = bndry_facenums[i]
        # get the local index of the vertices on the boundary face (local face number)
        # vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
        update_coords(mesh, bndry_i.element, real(vertices[itr][:,:,i]))
      end    # End for i = 1:nfaces
    end  # End for itr = 1:length(map.geom_faces)
  end # End

  commit_coords(mesh, sbp, opts)

  return nothing
end
=#

# complexify
import Base.isless

function isless(x::Number, y::Number)

  return isless(real(x), real(y))
end


end # End module FreeFormDeformation
