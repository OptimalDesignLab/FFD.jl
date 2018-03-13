module FFD


export AbstractMappingType, Mapping, PumiMapping
export PumiBoundingBox, calcKnot, controlPoint, calcParametricMappingLinear
export calcParametricMappingNonlinear, evalVolume, evalSurface
export writeControlPointsVTS, evaldXdControlPointProduct, evaldXdControlPointTransposeProduct

export initializeFFD, commitToPumi

export numLinearPlaneConstraints, countVarsLinearPlaneConstraints!, setLinearPlaneConstraints!
export numLinearCornerConstraints, countVarsLinearCornerConstraints!, setLinearCornerConstraints!
export numLinearStretchConstraints, countVarsLinearStretchConstraints!, setLinearStretchConstraints!
export numLinearRootConstraints, countVarsLinearRootConstraints!, setLinearRootConstraints!

export getControlPoints, setControlPoints

#push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))

using ArrayViews
using PumiInterface
using PdePumiInterface
using ODLCommonTools
using SummationByParts
using WriteVTK

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
  MPI.Init()
  atexit(finalizeMPI)
end



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

  This subtype of [`AbstractMapping`](@ref) enables the use of Free Form
  Deformation with a Pumi mesh.

  Users should *not* call this function directly, they should call
  [`initializeFFD`](@ref) instead.

  **Public Fields**

   * ndim: dimensionality of the mesh
   * cp_xyz: a 3 x `nctl[1]` x `nctl[2]` x `nctl[3]` array containing the
             x,y, and z coordinates of each control point.  Note that 
             the control point grid is always 3 dimensional, even when the
             mesh is 2 dimensional.  This array is read-only
    * n_face: a Pumi apf::Numbering object that numbers the face nodes
    * numFacePts: the number of points on the face
    * nctl: array of length 3, containing the number of control points in each
            direction.  For 2D meshes, the number of control points in the
            z direction must be 2.

  **Private Fields*

   * bc_nums: the array of boundary condition numbers that define the
              surface enveloped by the FFD control points
   * face_verts: vector of apf::MeshEntity* (Ptr{Void}) containing the
                 vertices in the order defined by `n_face`.  Length
                 `numFacePts`
   * map_cs: a PumiMapping{Complex128}, used for doing complex step when
             Tffd is real


  **Constructors**

   ***Constructor 1***
   
    This zero argument constructor produces an object with all fields
    zero.  Useful for constructing a dummy object

   ***Constructor 2***

    This constructor is the main constructor.  Its argument are:

     * ndims: number of dimensions (must be the same as mesh.dim
     * order: array of length 3 specifying the degree of the B-splines in
              each direction
     * nctl: array of length 3 specifying the number of control points in
             each direction
     * mesh: the Pumi mesh object
     * n_face: a Pumi apf::Numbering* hat numbers the surface points, default
               to C_NULL.  If unspecified, a new surface numbering will be
               created.
     * numFacePts: number of face points numbered by n_face.  Defaults to
                   zero.  This argument should only be specified in `n_face` is
                   as well.
     * bc_nums: array of integers specifying the boundary condition numbers
                (used in the construction of the `mesh`) that identify the
                surface to envelop in the FFD volume.  This is a keyword
                argument

  ***Constructor 3***

  Copy constructor.  Allows copying a PumiMapping.  The new one can have a
  different type parameter `Tffd`.  Arguments:

   * map: a PumiMapping
"""->

type PumiMapping{Tffd} <: AbstractMappingType{Tffd}

  ndim::Int                     # Mapping object to indicate 2D or 3D
  nctl::Array{Int, 1}   # Number of control points in each of the 3 dimensions
  order::Array{Int, 1}  # Order of B-spline in each direction

  xi::Array       # Paramaetric coordinates of input geometry
  cp_xyz::ROView{Tffd, 4, Array{Tffd, 4}}
  _cp_xyz::Array{Tffd, 4} # Cartesian coordinates of control points
  edge_knot::Array{Vector{Tffd}, 1}  # edge knot vectors
  bc_nums::Array{Int,1}
  cp_idx::Array{Int,4}  # index assigned to each CP coordinate

  # Working arrays
  aj::Array{Tffd, 3}
  dl::Array{Tffd, 2}
  dr::Array{Tffd, 2}
  work::Array{Tffd, 4}

  n_face::Ptr{Void}  # apf::Numbering*, numbering the face points
  numFacePts::Int  # number of face points
  face_verts::Array{Ptr{Void}, 1}  # face points, in order

  map_cs::PumiMapping{Complex128}  # map for complex stepping even when Tffd
                                   # is real
  vertices::Array{Tffd, 2}  # dim x numFacePts array, used by map_cs

  function PumiMapping(ndim::Int, order::AbstractArray{Int,1},
                       nctl::AbstractArray{Int,1}, mesh::AbstractMesh,
                       n_face::Ptr{Void}=C_NULL, numFacePts=0;
                       bc_nums::AbstractArray{Int,1}=[0])

    #TODO: remove dims, use mesh.dim instead
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
    map.order = order
    map.nctl = nctl
    map._cp_xyz = zeros(3, nctl[1], nctl[2], nctl[3])
    map.cp_xyz = ROView(map._cp_xyz)

    map.edge_knot = Array(Vector{Tffd}, 3)
    for i = 1:3
      map.edge_knot[i] = zeros(Tffd, nctl[i]+order[i])
    end

    # if the user didn't provide a surface point numbering, make one ourselves
    if n_face == C_NULL
      numFacePts, n_face, face_verts = numberSurfacePoints(mesh, bc_nums)
    else  # get the face_verts array
      face_verts = Array(Ptr{Void}, numFacePts)
      for vert in mesh.verts
        n_v = getNumberJ(n_face, vert, 0, 0)
        if n_v <= numFacePts
          face_verts[n_v] = vert
        end
      end
    end
    map.n_face = n_face
    map.numFacePts = numFacePts
    map.face_verts = face_verts
    map.xi = zeros(Tffd, 3, numFacePts)
    map.bc_nums = bc_nums


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

    # do this in initializeFFD because map isn't fully initialized yet
#    map.map_cs = PumiMapping{Complex128}(map)
    map.vertices = Array(Tffd, 0, 0)

    return map
  end

  function PumiMapping()
    map = new()
    map.ndim = 0
    map.nctl = Int[]
    map.order = Int[]
    map.xi = zeros(Tffd, 0, 0)
    map._cp_xyz = zeros(Tffd, 0, 0, 0, 0)
    map.cp_xyz = ROView(map._cp_xyz)
    map.edge_knot = Array(Vector{Tffd}, 0)
    map.bc_nums = Int[]
    map.cp_idx = zeros(Tffd, 0, 0, 0, 0)
    map.aj = zeros(Tffd, 0, 0, 0)
    map.dl = zeros(Tffd, 0, 0)
    map.dr = zeros(Tffd, 0, 0)
    map.work = zeros(Tffd, 0, 0, 0, 0)

    map.n_face = C_NULL
    map.numFacePts = 0
    map.face_verts = Ptr{Void}[]
    map.map_cs = PumiMapping{Tffd}(map)
    map.vertices = Array(Tffd, 0, 0)

    return map
  end

  # constructs a new mapping from a fully initialized one, with possibly
  # different parameter Tffd
  function PumiMapping(map_old::PumiMapping)

    map = new()
    map.ndim = map_old.ndim
    map.nctl = copy(map_old.nctl)
    map.order = copy(map_old.order)
    map.xi = zeros(Tffd, size(map_old.xi)); copy!(map.xi, map_old.xi)
    map._cp_xyz = zeros(Tffd, size(map_old.cp_xyz)); copy!(map._cp_xyz, map_old._cp_xyz)
    map.cp_xyz = ROView(map._cp_xyz)
    map.edge_knot = Array(Vector{Tffd}, length(map_old.edge_knot))
    for i=1:length(map.edge_knot)
      map.edge_knot[i] = copy(map_old.edge_knot[i])
    end
    map.bc_nums = copy(map_old.bc_nums)
    map.cp_idx = copy(map_old.cp_idx)

    map.aj = zeros(Tffd, size(map_old.aj)); copy!(map.aj, map_old.aj)
    map.dl = zeros(Tffd, size(map_old.dl)); copy!(map.dl, map_old.dl)
    map.dr = zeros(Tffd, size(map_old.dr)); copy!(map.dr, map_old.dr)
    map.work = zeros(Tffd, size(map_old.work)); copy!(map.work, map_old.work)

    map.n_face = map_old.n_face
    map.numFacePts = map_old.numFacePts
    map.face_verts = map_old.face_verts  # reference, don't copy for this one
                                         # because it is large
    # leave map_cs uninitialized

    map.vertices = Array(Tffd, map_old.ndim, map.numFacePts)
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
* `bc_nums`: (Optional Argument) If partial geometry embedded, give the
                face/edge numbers of the embedded gemetry
* `n_face`: an apf::Numbering* (Ptr{Void}) numbering the face nodes, if
            not supplied a new numbering will be created
* `numFacePts`: the number of face points on the surface defined by
                `bc_nums` and numbered by `n_face`

"""->

function initializeFFD{Tmsh}(mesh::AbstractMesh{Tmsh}, sbp::AbstractSBP,
                       order::AbstractArray{Int,1},
                       nControlPts::AbstractArray{Int,1},
                       offset::AbstractArray{Float64,1},
                       bc_nums::AbstractArray{Int,1}=[0],
                       n_face::Ptr{Void}=C_NULL,
                       numFacePts::Integer=0)

  # Create Mapping object
  ndim = mesh.dim
  map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, n_face, numFacePts,
                              bc_nums=bc_nums)

  # Create knot vector
  calcKnot(map)

  # Create Bounding box
  ffd_box = PumiBoundingBox{Tmsh}(map, mesh, sbp, offset)

  # Control points
  controlPoint(map, ffd_box)

  # Populate map.xi
  calcParametricMappingNonlinear(map, ffd_box, mesh, bc_nums)

  # now that map is fully initialized, copy it for the complex step calculations
  map.map_cs = PumiMapping{Complex128}(map)
  return map, ffd_box
end

"""
  This function copies the supplied array into map._cp_xyz, without checking
  the dimensionality first.

  **Inputs**

   * map: a PumiMapping object
   * cp_xyz: array, same size as map.cp_xyz, with new control point coordinates.
"""
function _setControlPoints{T}(map::PumiMapping, cp_xyz::AbstractArray{T, 4})

  for i=1:4
    @assert size(cp_xyz, i) == size(map._cp_xyz, i)
  end

  for i=1:length(cp_xyz)
    map._cp_xyz[i] = cp_xyz[i]
  end

  return nothing
end

"""
  This function updates the control point locations.  Users should not modify
  map.cp_xyz directly, they should use this function instead.  Derivative
  values will be incorrect if map.cp_xyz is modified.

  **Inputs**

   * map: a PumiMapping
   * cp_xyz: array with new control point locations.  In 2D, this array 
             should be 2 x map.nctl[1] x map.nctl[2].  In 3D it should be
             3 x map.nctl[1] x map.nctl[2] x map.nctl[3]

  **Implementaton Notes**

  The internals of FFD always assume a 3 dimensional mesh, but for running
  2D problems having 2 control points along the z axis doesn't make sense.
  In 2D, this function updates both planes of control points simultaneously.
"""
function setControlPoints{T}(map::PumiMapping, cp_xyz::AbstractArray{T, 4})

  @assert map.ndim == 3

  # copy everything over
  _setControlPoints(map, cp_xyz)

  return nothing
end

function setControlPoints{T}(map::PumiMapping, cp_xyz::AbstractArray{T, 3})

  @assert map.ndim == 2
  @assert size(cp_xyz, 1) == 2
  @assert size(cp_xyz, 2) == size(map.cp_xyz, 2)
  @assert size(cp_xyz, 3) == size(map.cp_xyz, 3)
  @assert size(map.cp_xyz, 4) == 2

  # set cp_xyz for both z = zmin and z = zmax
  # note that the derivative calculations need to take this into account

  for i=1:size(cp_xyz, 3)
    for j=1:size(cp_xyz, 2)
      for k=1:size(cp_xyz, 1)
        map._cp_xyz[k, j, i, 1] = cp_xyz[k, j, i]
        map._cp_xyz[k, j, i, 2] = cp_xyz[k, j, i]
      end
    end
  end

  return nothing
end

"""
  This function returns a copy of the current control point locations, suitable
  to be passed into [`setControlPoints`](@ref)  See that function for details

  Note that this function is type unstable.  Users are generally better off
  allocating their own array.

  **Inputs**

   * map: PumiMapping

  **Outputs**

   * array of control point locations
"""
function getControlPoints(map::PumiMapping)

  if map.ndim == 2
    return map._cp_xyz[1:2, :, :, 1]
  else
    return copy(map._cp_xyz)
  end

  return nothing
end

"""
  This function checks if an array formatted similarly to map.cp_xyz
  has the same coordinates for the z = zmin and z = zmax.  This is the
  condition enforced by [`setControlPoints`](@ref) for 2D arrays.

  **Inputs**

   * cp_xyz: array like cp_xyz
"""
function check2DSymmetry{T}(cp_xyz::AbstractArray{T, 4}, tol=1e-13)

  @assert size(cp_xyz, 4) == 2

  nonsym = false
  for i=1:size(cp_xyz, 3)
    for j=1:size(cp_xyz, 2)
      for k=1:2
        flag = abs(cp_xyz[k, j, i, 1] - cp_xyz[k, j, i, 2]) > tol
        nonsym = nonsym || flag
      end
    end
  end

  @assert !nonsym "cp_xyz is non-symmetric about the x-y plane.  Did you modify PumiMapping.cp_xyz directly?"

  return nothing
end


#=
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
=#

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


end # module
