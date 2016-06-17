# startup.jl
#
# Startup script for FFD. This will change as the code progresses

# Includes
include("mapping.jl")

using ArrayViews

ndim = 3
order = [2,2,2]  # Order of B-splines in the 3 directions
nControlPts = [3,3,3]
nnodes = [3,3,3]  # Number of nodes of the FE grid that need to be mapped

map = LinearMapping(ndim, order, nControlPts, nnodes)
calcKnot(map)  # Create knot vectors

for i = 1:3
  println("map.edge_knot[$i] = ", map.edge_knot[i])
end

# Create Bounding Box
offset = [0.5,0.5,0.5]  # offset for the bounding box
geom_bounds = [1. 1. 1.;3. 3. 3.]
box = BoundingBox(ndim, geom_bounds, offset)
println("origin = $(box.origin)")

# Unit vectors
# box.unitVector = eye(Float64, 3)

println("box.lengths = $(box.lengths)")
println("box.box_bound = $(box.box_bound)\n\n")

#------------------------------------
# Representative obeject
# nnodes = [5,3,4]
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
#=
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_xyz[$i,$j,$k,:] = ", nodes_xyz[i,j,k,:])
    end
    println('\n')
  end
end
=#
#----------------------------------------

calcParametricMappingLinear(map, box, nodes_xyz)
#=
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_stu[$i,$j,$k,:] = ", map.xi[i,j,k,:])
    end
  end
end
=#
# Define the control points
controlPoint(map,box)

for k = 1:map.nctl[3]
  for j = 1:map.nctl[2]
    for i = 1:map.nctl[1]
      println("cp_xyz[$i,$j,$k,:] = ", round(map.cp_xyz[i,j,k,:],2) )
    end
    println('\n')
  end
end

Vol = zeros(nodes_xyz)
# evalVolume(map, Vol)
#=
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_xyz_err[$i,$j,$k,:] = ", nodes_xyz[i,j,k,:]-Vol[i,j,k,:])
    end
    println('\n')
  end
  println('\n')
end
=#

# Check linear scaling (Translation)
#=
fill!(Vol, 0.0)
# move the x coordinates of all the control by 4 units
map.cp_xyz[:,:,:,3] += 2
evalVolume(map, Vol)
# Output the values to check for errors
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println( "Vol_err[$i,$j,$k,:] = ", round(Vol[i,j,k,:]-nodes_xyz[i,j,k,:],2) )
    end
    println('\n')
  end
  println('\n')
end
=#

# Check angular rotation
# Rotate control points by 90 degrees about Z axis
gamma = 0.5*pi # angle in radians by which ctls are rotated
rotx = [1. 0. 0.;0. cos(gamma) -sin(gamma);0. sin(gamma) cos(gamma)]
roty = [cos(gamma) 0. sin(gamma);0. 1. 0.;-sin(gamma) 0. cos(gamma)]
rotz = [cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0;0 0 1]
controlPoint(map,box) # Recreate unperturbed control points
# Perform the rotation
for k = 1:map.nctl[3]
  for j = 1:map.nctl[2]
    for i = 1:map.nctl[1]
      map.cp_xyz[i,j,k,:] = rotx*map.cp_xyz[i,j,k,:]
    end
  end
end


println("Rotated coordinates\n")#=
for k = 1:map.nctl[3]
  for j = 1:map.nctl[2]
    for i = 1:map.nctl[1]
      println("cp_xyz[$i,$j,$k,:] = ", round(map.cp_xyz[i,j,k,:],2) )
    end
    println('\n')
  end
end=#

evalVolume(map, Vol)
# Check Output
#=
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println( "Vol[$i,$j,$k,:] = ", round(Vol[i,j,k,:],2) )
    end
    println('\n')
  end
  println('\n')
end=#
