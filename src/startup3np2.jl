# startup3.jl
# Merging mesh movement and FreeFormDeformation

push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

gc()

using PdePumiInterface
using SummationByParts
using ODLCommonTools
using ArrayViews
using Utils
using MPI
using MeshMovement
using FreeFormDeformation

# MPI Declarations
MPI.Init()
comm = MPI.COMM_WORLD
comm_world = MPI.MPI_COMM_WORLD
comm_self = MPI.COMM_SELF
my_rank = MPI.Comm_rank(comm)
comm_size = MPI.Comm_size(comm)

opts = PdePumiInterface.get_defaults()
# 2D mesh
opts["order"] = 1
opts["dimensions"] = 2
opts["use_DG"] = true
opts["operator_type"] = "SBPOmega"
opts["dmg_name"] = "./mesh_files/2DAirfoil_chord1.dmg"
opts["smb_name"] = "./mesh_files/2DAirfoil_chord1.smb"
opts["numBC"] = 2

# For 2DAirfoil
opts["BC1"] = [8]
opts["BC1_name"] = "FarField"
opts["BC2"] = [5]
opts["BC2_name"] = "Airfoil"

opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
opts["jac_type"] = 2

sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

geom_faces = opts["BC2"]

# Free Form deformation parameters
ndim = 2
order = [4,4,2]  # Order of B-splines in the 3 directions
nControlPts = [6,6,2]

# Create Mapping object
ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false, geom_faces=geom_faces)

# Create knot vector
calcKnot(ffd_map)

# Create Bounding box
offset = [0., 0., 0.5] # No offset in the X & Y direction
ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)

# Control points
controlPoint(ffd_map, ffd_box)

# Populate map.xi
calcParametricMappingNonlinear(ffd_map, ffd_box, mesh, geom_faces)

# Prep MeshWarping
volNodes = zeros(Tmsh, 3, mesh.numVert)
MPI.Barrier(comm)
println("size mesh.vert_coords = $(size(mesh.vert_coords))")
for i = 1:mesh.numEl
  for j = 1:size(mesh.vert_coords,2)
    # Get the vertex numbering on the portion of mesh owned by the processor
    local_vertnum = mesh.element_vertnums[j,i]
    volNodes[1:2,local_vertnum] = mesh.vert_coords[:,j,i]# mesh.element_vertnums
  end
end
orig_volNodes = deepcopy(volNodes)

# Get the bool array of all the surface vertex coordinates
surfaceVtx = zeros(Int32, mesh.numVert)
for (bindex, bndry) in enumerate(mesh.bndryfaces)
  vtx_arr = mesh.topo.face_verts[:, bndry.face]
  for j = 1:length(vtx_arr)
    local_vertnum = mesh.element_vertnums[vtx_arr[j], bndry.element]
    surfaceVtx[local_vertnum] = 1
  end
end

# Prepare the wall coordinates array for mesh warping
nwall_faces = zeros(Int,length(geom_faces))
vtx_per_face = mesh.dim # only true for simplex elements
for itr = 1:length(geom_faces)
  geom_face_number = geom_faces[itr]
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
  nwall_faces[itr] = length(bndry_facenums)
end      # End for itr = 1:length(geomfaces)

# Populate wallCoords
wallCoords = zeros(3, sum(nwall_faces)*vtx_per_face)
ctr = 1 # Counter for wall coordinates
for itr = 1:length(geom_faces)
  geom_face_number = geom_faces[itr]
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
    # get the local index of the vertices
    vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
    for j = 1:length(vtx_arr)
      wallCoords[1:2,ctr] = mesh.vert_coords[:,vtx_arr[j],bndry_i.element]
      ctr += 1
    end  # End for j = 1:length(vtx_arr)
  end    # End for i = 1:nfaces
end

# Warp parameters
param = warpParam(nLocalVolNodes=mesh.numVert)

# MPI Information for warping
mpiVar = MPIInfo(warp_comm_world=comm_world, myID=my_rank, nProc=comm_size)

# Populate entries of param
param.aExp = 2.
param.bExp = 2.
param.LdefFact = 1
param.alpha = 0.2
param.symmTol = 1e-4
param.errTol = 1e-4
param.cornerAngle = 30.
param.zeroCornerRotations = false
param.useRotations = false
param.evalMode = 5

# Symmetry Information
# A symmetry plane is defined a point and normal. It is necessary to ensure that
# Each eplane is independent of the other. Look at checkPlane by Gaetan in
# UnstructuredMesh.Py. The way this is done in the data structure is that, every
# plane is defined by a (3,2) array. [:,1] defines the point and [:,2] defines
# the unit normal vector.
symmetryPlanes = zeros(Float64,3,2)
flatWarpSurfPts = reshape(wallCoords, 3*size(wallCoords,2))
faceSizes = mesh.dim*ones(Int32,sum(nwall_faces))

initializeParameters(param, mpiVar, symmetryPlanes)

initializeWarping(param, volNodes, wallCoords, faceSizes, comm)

# New surface coordinates
#=
# Camber the shit
ffd_map.cp_xyz[2,1:2,1:4,:] -= 0.01*4
ffd_map.cp_xyz[2,7,1:4,:] -= 0.015*4
ffd_map.cp_xyz[2,4:6,5:7,:] += 0.01*4
=#

# Translate control points along x & y by  5 units
ffd_map.cp_xyz[1,:,:,:] += 0.1*4
ffd_map.cp_xyz[2,:,:,:] += 0.5*4


# Rotation matrix
theta = -20*pi/180  # Rotate wall coordinates by 10 degrees
rotMat = [cos(theta) -sin(theta) 0
          sin(theta) cos(theta)  0
          0          0           1] # Rotation matrix

# Rotate the control points
for k = 1:ffd_map.nctl[3]
  for j = 1:ffd_map.nctl[2]
    for i = 1:ffd_map.nctl[1]
      ffd_map.cp_xyz[:,i,j,k] = rotMat*ffd_map.cp_xyz[:,i,j,k]
    end
  end
end

evalSurface(ffd_map, mesh)

# Update wallCoords
ctr = 1 # Counter for wall coordinates
for itr = 1:length(geom_faces)
  geom_face_number = geom_faces[itr]
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
    # get the local index of the vertices
    vtx_arr = mesh.topo.face_verts[:,bndry_i.face]
    for j = 1:length(vtx_arr)
      wallCoords[1:2,ctr] = mesh.vert_coords[:,vtx_arr[j],bndry_i.element]
      #println("wallCoords[:, $ctr] = $(wallCoords[:,ctr])")
      ctr += 1
    end  # End for j = 1:length(vtx_arr)
  end    # End for i = 1:nfaces
end

MPI.Barrier(comm)
warpMesh(param, volNodes, wallCoords)
MPI.Barrier(comm)
for i = 1:mesh.numEl
  for j = 1:size(mesh.vert_coords,2)
    # Get the vertex numbering on the portion of mesh owned by the processor
    local_vertnum = mesh.element_vertnums[j,i]
    mesh.vert_coords[:,j,i] = volNodes[1:2,local_vertnum] # mesh.element_vertnums
  end
end

for i = 1:mesh.numEl
  update_coords(mesh, i, mesh.vert_coords[:,:,i])
end
commit_coords(mesh, sbp)
writeVisFiles(mesh, "TestNewAPI")
MPI.Finalize()
