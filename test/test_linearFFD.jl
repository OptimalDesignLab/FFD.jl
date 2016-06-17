# test components
include("mapping.jl")
include("knot.jl")
include("bounding_box.jl")
include("linear_mapping.jl")
include("control_point.jl")
include("span.jl")
include("b-splines.jl")
include("evaluations.jl")

using ArrayViews
using FactCheck

ndim = 3
order = [4,4,4]  # Order of B-splines in the 3 directions
nControlPts = [4,4,4]
nnodes = [5,3,4]  # Number of nodes of the FE grid that need to be mapped
map = Mapping(ndim, order, nControlPts, nnodes) # Create Mapping Object

# Create BoundingBox object
offset = [0.,0.,0.]  # offset for the bounding box
geom_bounds = [1. 1. 1.;5. 3. 4.]
box = BoundingBox(ndim, geom_bounds, offset)

facts("--- Checking Mapping object ---") do

  @fact map.ndim --> 3
  @fact map.nctl --> [4,4,4]
  @fact map.order --> [4,4,4]
  @fact map.numnodes --> [5,3,4]

  context("--- Checking Knot calculations ---") do

    calcKnot(map)
    for i = 1:map.ndim
      @fact map.edge_knot[i] --> roughly([0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0],
                                         atol=1e-15)
    end

  end  # End context("--- Checking Knot calculations ---")

end  # End facts("--- Checking Mapping object ---") do

facts("--- Checking BoundingBox ---") do

  @fact box.ndim --> 3
  @fact box.origin --> [1., 1., 1.]
  @fact box.unitVector --> [1. 0. 0.;0. 1. 0.;0. 0. 1.]
  @fact geom_coord --> [1. 1. 1.;5. 3. 4.]
  @fact box.offset --> [0.,0.,0.]

end  # End facts("--- Checking BoundingBox ---") do

# Create Test mesh
nodes_xyz = zeros(nnodes[1], nnodes[2], nnodes[3], 3)
origin = [1,1,1]
incz = 0.0
for k = 1:nnodes[3]
  incy = 0.0
  for j = 1:nnodes[2]
    incx = 0.0
    for i = 1:nnodes[1]
      nodes_xyz[i,j,k,1] = origin[1] + incx
      nodes_xyz[i,j,k,2] = origin[2] + incy
      nodes_xyz[i,j,k,3] = origin[3] + incz
      incx += 1
    end
    incy += 1
  end
  incz += 1
end

facts("--- Checking Linear Mapping ---") do


end # End facts("--- Checking Linear Mapping ---") do
