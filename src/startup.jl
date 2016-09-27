# startup.jl
#
# Startup script for FFD. This will change as the code progresses

# Includes
push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))

using PdePumiInterface
using SummationByParts
using ODLCommonTools
using ArrayViews
using MPI
using MeshMovement

using FreeFormDeformation

opts = PdePumiInterface.get_defaults()
# 2D mesh
opts["numBC"] = 2
opts["BC1"] = [8,11,14,17]
opts["BC1_name"] = "FarField"
opts["BC2"] = [5]
opts["BC2_name"] = "Airfoil"
opts["coloring_distance"] = 0 # For CG Mesh

Tmsh = Float64
dmg_name = ".null"
smb_name = "./mesh_files/gvortex1.smb"
order = 1
dofpernode = 1

# SBP & Mesh Parameters
# For SBP Gamma
reorder = true
internal = false
shape_type = 1

# create linear sbp operator
sbp = TriSBP{Float64}(degree=order, reorder=reorder, internal=internal)
sbpface = TriFace{Float64}(order, sbp.cub, sbp.vtx)
# create linear mesh with 4 dof per node

println("constructing CG mesh")
mesh = PumiMesh2{Tmsh}(dmg_name, smb_name, order, sbp, opts, sbpface;
       dofpernode=dofpernode, coloring_distance=opts["coloring_distance"],
       shape_type=shape_type)
#
#
# Free Form deformation parameters
ndim = 2
order = [4,4,2]  # Order of B-splines in the 3 directions
nControlPts = [4,4,2]
mesh_info = Int[sbp.numnodes, mesh.numEl]

ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh_info)
calcKnot(ffd_map)
println("ffd_map.edge_knot = \n", ffd_map.edge_knot)

# Create Bounding box
offset = [0., 0., 0.5]
ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)

# println("ffd_box.origin = ", ffd_box.origin)


controlPoint(ffd_map, ffd_box)
#=
for k = 1:ffd_map.nctl[3]
  for j = 1:ffd_map.nctl[2]
    for i = 1:ffd_map.nctl[1]
      println("cp_xyz[3,$i,$j,$k] = ", ffd_map.cp_xyz[3,i,j,k])
    end
    println('\n')
  end
end
=#
# Create Linear Mapping
calcParametricMappingNonlinear(ffd_map, ffd_box, mesh)
#println("ffd_map.xi = \n", ffd_map.xi)
#=
# Translate control points along x & y by  5 units
ffd_map.cp_xyz[1,:,:,:] += 2
ffd_map.cp_xyz[2,:,:,:] += 3
=#
# Check angular rotation
# Rotate control points by 90 degrees about Z axis
gamma = 0.5*pi # angle in radians by which ctls are rotated
rotx = [1. 0. 0.;0. cos(gamma) -sin(gamma);0. sin(gamma) cos(gamma)]
roty = [cos(gamma) 0. sin(gamma);0. 1. 0.;-sin(gamma) 0. cos(gamma)]
rotz = [cos(gamma) -sin(gamma) 0;sin(gamma) cos(gamma) 0;0 0 1]
# Perform the rotation
for k = 1:ffd_map.nctl[3]
  for j = 1:ffd_map.nctl[2]
    for i = 1:ffd_map.nctl[1]
      ffd_map.cp_xyz[:,i,j,k] = rotz*ffd_map.cp_xyz[:,i,j,k]
    end
  end
end

evalVolume(ffd_map, mesh)

# println("mesh.coords = \n", mesh.coords)

for i = 1:mesh.numEl
  update_coords(mesh, i, mesh.coords[:,:,i])
end
PumiInterface.writeVtkFiles("Translation", mesh.m_ptr)


# geom_bounds = zeros(2,3)
# FreeFormDeformation.calcGeomBounds(mesh.coords, geom_bounds)

#=
map = Mapping(ndim, order, nControlPts, nnodes)
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
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_xyz[$i,$j,$k,:] = ", nodes_xyz[i,j,k,:])
    end
    println('\n')
  end
end
#----------------------------------------

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

# calcParametricMappingLinear(map, box, nodes_xyz)
calcParametricMappingNonlinear(map, box, nodes_xyz)

for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_stu[$i,$j,$k,:] = ", map.xi[i,j,k,:])
    end
  end
end


Vol = zeros(nodes_xyz)
evalVolume(map, Vol)

for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println("nodes_xyz_err[$i,$j,$k,:] = ", round(Vol[i,j,k,:]-nodes_xyz[i,j,k,:],2))
    end
    println('\n')
  end
  println('\n')
end


# Check linear scaling (Translation)
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


println("Rotated coordinates\n")
for k = 1:map.nctl[3]
  for j = 1:map.nctl[2]
    for i = 1:map.nctl[1]
      println("cp_xyz[$i,$j,$k,:] = ", round(map.cp_xyz[i,j,k,:],2) )
    end
    println('\n')
  end
end

evalVolume(map, Vol)
# Check Output
for k = 1:nnodes[3]
  for j = 1:nnodes[2]
    for i = 1:nnodes[1]
      println( "Vol[$i,$j,$k,:] = ", round(Vol[i,j,k,:],2) )
    end
    println('\n')
  end
  println('\n')
end
=#
