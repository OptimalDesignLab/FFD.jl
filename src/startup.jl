# startup.jl
#
# Startup script for FFD. This will change as the code progresses

# Includes
include("mapping.jl")
include("knot.jl")
include("bounding_box.jl")
# include("linear_mapping.jl")

ndim = 3
order = [4,4,4]  # Order of B-splines in the 3 directions
nControlPts = [8,8,8]
nnodes = [11,11,11]  # Number of nodes of the FE grid that need to be mapped

map = Mapping(ndim, order, nControlPts, nnodes)

# Create Bounding Box
offset = [0.,0.,0.]  # offset for the bounding box
geom_bounds = [3 3 3;10 6 10]
box = BoundingBox(ndim, geom_bounds, offset)
println("origin = $(box.origin)")

# Unit vectors
box.unitVector = eye(Float64, 3)

x = [10.,4.5,10.]
pX = zeros(3)
println("box.lengths = $(box.lengths)")
println("box.box_bound = $(box.box_bound)\n\n")

#------------------------------------
# Representative obeject
nnodes = [4,2,3]
nodes_xyz = zeros(nnodes[1], nnodes[2], nnodes[3], 3)
origin = [1,1,1]
incx = 0.0
incy = 0.0
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

# Output the values to check for errors
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_xyz[$i,$j,$k,:] = ", nodes_xyz[i,j,k,:])
    end
    println('\n')
  end
  println('\n')
end


#----------------------------------------
#=
linearMap(map, box, x, pX)

println("pX = $pX")
=#
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
