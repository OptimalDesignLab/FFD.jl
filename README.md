# FFD.jl
A Julia repository containing code for Free Form Deformation using B-splines.
You only need to clone the repository to use it. There is no need to build
anything. Please read the `README.md` file to understand how to use the
repository. Within each of the files, there is detailed documentation of each
and every function.

The repository has the following structure

* **`FFD.jl`**
  * **`src`**
    * `bounding_box.jl`
    * `b-splines.jl`
    * `control_point.jl`
    * `evaluations.jl`
    * `extra.jl`
    * `knot.jl`
    * `mapping_functions.jl`
    * `mapping.jl`
    * `param_functions.jl`
    * `plot.jl`
    * `span.jl`
    * `startup.j`l
  * **`test`**
    * `runtest.jl`
    * `test_b-splines.jl`
    * `test_linearFFD.jl`
  * `README.md`

The code was developed from the following references

1.  Les Piegl and Wayne Tiller. 1997. *The NURBS Book (2nd Ed.)*.
Springer-Verlag New York, Inc., New York, NY, USA.
2.  Carl de Boor. 1978. *A Prcatical Guide to Splines*. Springer-Verlag
New York, Inc., New York, NY, USA.
3.  Thomas W. Sederberg and Scott R. Parry. 1986. Free-form deformation of solid
geometric models. In *Proceedings of the 13th annual conference on Computer
graphics and interactive techniques* (SIGGRAPH '86), David C. Evans and Russell
J. Athay (Eds.). ACM, New York, NY, USA, 151-160.
DOI=http://dx.doi.org/10.1145/15922.15903
4. Prof. Hicken's Mapping_Mod.f90
5. University of Michigan Free-form Deformation Toolbox

The file `startup.jl` provides a sample method to call the functions within the
repository. There are two major composite types within the repository that needed
for FFD, viz. `Mapping` and `BoundingBox`.

`Mapping` type essentially constructs a mapping object which not only creates the
mapping between the physical and parametric coordinate systems, but also store
the knot vectors, control points and other variables that are necessary for
construction of splines and FFD in general. Those details can be found in the
documentation. An important thing to note here is the notation of the physical
and parametric spaces. While the physical spaces is denoted by *(x,y,z)*, the
parametric spaces are represented as follows:

* within the functions that compute B-spline basis functions, knot span index,
B-spline curve points and B-spline volume points, the parametric coordinates are
represented as *(u,v,w)* taken from Ref[1].
* The parametric coordinates are represented as *(ξ,η,ζ)* within the Mapping
object. The notation is taken from Ref[4].
* Within the functions that compute the linear and non-linear mapping, the
parametric coordinates are referred to as *(s,t,u)*

The intention of doing this presently is to ensure easy comprehension of the
code if modifications needed referencing the literature. They may be changed in
the future.

`BoundingBox` type contains certain geometric information about the FFD volume
in which the a geometry is embedded. Its capabilities are restricted presently
and will grow as complicated FFD volumes are considered in the future.

Two types of mapping can be constructed between the physical and the parametric
space, viz. *linear mapping* and *nonlinear mapping*. A linear mapping will
suffice if the FFD volume considered is a regular parallelepiped. However, if
the FFD volume is an irregular shape, a nonlinear mapping is required which
performs a Newton's solve to obtain the parametric coordinates.

A sample step by step procedure for using FFD would be
```
# Import geometry

# Create Mapping Object
ndim = 3
order = [2,2,2]  # Order of B-splines in the 3 directions
nControlPts = [3,3,3]
map = Mapping(ndim, order, nControlPts, nnodes)

# Create BoundingBox object
offset = [0.5,0.5,0.5]  # offset for the bounding box
geom_bounds = [1. 1. 1.;3. 3. 3.]
box = BoundingBox(ndim, geom_bounds, offset)

calcKnot(map)         # Create knot vectors
controlPoint(map,box) # Create Control Points for FFD volume

# Create Linear or nonlinear mapping
if type_of_map == "linear"
  calcParametricMappingLinear(map, box, embedded_geometry_node_coordinates)
else
  calcParametricMappingNonlinear(map, box, embedded_geometry_node_coordinates)
end

# Array for storing the computed embedded geometry coordinates within the FFD volume
FFD_vol = zeros(embedded_geometry_node_coordinates)

# Manipulate Control points

# Function to evaluate the computed embedded geometry coordinates within the FFD volume
evalVolume(map, Vol)

## VOILA!!!
```

When using FFD with a Pumi mesh, an example usage is
```
  mesh, sbp, opts = getTestMesh()  # get mesh, sbp, options dictionary from somewhere
  Tmsh = eltype(mesh.jac)  # get the element type of the mesh arrays

  # create the FFD around the specified surface
  order = [3,3,3]  # Order of B-splines in the 3 directions
  nControlPts = [6,6,6]
  offset = [0.25, 0.25, 0.25]
  bc_nums = [1]  # envelop the surface that boundary condition 1 is applied
                 # to in the FFD volume

  map, box = initializeFFD(mesh, sbp, order, nControlPts, offset, 
                           bc_nums)


  # evaluate surface point locations
  vertices_orig = zeros(Tmsh, mesh.dim, map.numFacePts)
  evalSurface(map, mesh, vertices_orig)

  # change control points locations
  map.cp_xyz[1,:,:,:] += 0.1
  map.cp_xyz[2,:,:,:] += 0.2
  map.cp_xyz[3,:,:,:] += 0.3

  # get updated coordinates
  vertices_new = zeros(Tmsh, mesh.dim, map.numFacePts)
  evalSurface(map, mesh, vertices_new)

  # compute the Jacobian-vector product with a random vector
  Xcp_dot = rand(map.cp_xyz)
  Xs_dot = zeros(vertices_new)
  evaldXdControlPointProduct(map, mesh, Xcp_dot, Xs_dot)
  # Xs_dot now contains the result

  # evaluate the transposed Jacobian-vector product with a random vector
  Xcp_bar = zeros(map.cp_xyz)
  Xs_bar = rand(vertices_new)
  evaldXdControlPointTransposeProduct(map, mesh, Xs_bar, Xcp_bar)
  # Xcp_bar now has contains the results
```

See the docstrings of these functions for more details.

[![Build Status](https://travis-ci.org/OptimalDesignLab/FFD.jl.svg?branch=master)](https://travis-ci.org/OptimalDesignLab/FFD.jl)
