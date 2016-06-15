# startup.jl
#
# Startup script for FFD. This will change as the code progresses

# Includes
include("mapping.jl")
include("knot.jl")
include("bounding_box.jl")
include("linear_mapping.jl")

ndim = 3
order = [4,4,4]  # Order of B-splines in the 3 directions
nControlPts = [8,8,8]
nnodes = [11,11,11]  # Number of nodes of the FE grid that need to be mapped

map = Mapping(ndim, order, nControlPts, nnodes)

# Create Bounding Box
offset = [0.,0.,0.]  # offset for the bounding box
geom_bounds = [5 5 5;10 10 10]
box = BoundingBox(ndim, geom_bounds, offset)
println("origin = $(box.origin)")

# Unit vectors
box.unitVector = 5*eye(Float64, 3)

x = [10.,10.,10.]
pX = zeros(3)

linearMap(map, box, x, pX)

println("pX = $pX")

# Populate S, T, U direction vectors for the

#=
println(map.ndim)
println(map.nctl)

#  knot vectors
calcKnot(map)
# println(map.edge_knot)

U = zeros(3,3,3)
V = zeros(U)
W = zeros(U)

for i = 1:3
  for j = 1:3
    for k = 1:3
      U[i,j,k] = 0.5*i
      V[i,j,k] = 0.5*j
      W[i,j,k] = 0.5*k
    end
  end
end

for i = 1:3
  println("\nU[:,:,$i] = \n", U[:,:,i])
end
println("\nV = \n", V)
println("\nW = \n", W)
=#
