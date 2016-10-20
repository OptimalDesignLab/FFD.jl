# startup2.jl
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
opts["dmg_name"] = "./mesh_files/2D_Airfoil.dmg"
opts["smb_name"] = "./mesh_files/2D_Airfoil.smb"
opts["numBC"] = 2
opts["BC1"] = [8,11,14,17]
opts["BC1_name"] = "FarField"
opts["BC2"] = [5]
opts["BC2_name"] = "Airfoil"
opts["coloring_distance"] = 2 # 0 For CG Mesh 2 for DG Mesh
opts["jac_type"] = 2

sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

geom_faces = opts["BC2"]
println("geom_faces = $geom_faces")

# Free Form deformation parameters
ndim = 2
order = [4,4,2]  # Order of B-splines in the 3 directions
nControlPts = [4,4,2]

# Copy the original vertex coordinates
orig_vert_coords = deepcopy(mesh.vert_coords)


ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false,
geom_faces=[5])
println("ffd_map.geom_faces = $(ffd_map.geom_faces)")
println("ffd_map.xi[1] = $(size(ffd_map.xi[1]))")
println("size_mesh.vert_coords = $(size(mesh.vert_coords))")
calcKnot(ffd_map)
println("knot vectors = \n", ffd_map.edge_knot)

# Create Bounding box
offset = [0., 0., 0.5] # No offset in the X & Y direction
ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)
println("box origin = $(ffd_box.origin)")
println("box bounds = \n", ffd_box.box_bounds)

# Control points
controlPoint(ffd_map, ffd_box)

# Populate map.xi
calcParametricMappingLinear(ffd_map, ffd_box, mesh, geom_faces)

#=
# Check if xi is being computed correctly
err_ctr = 0
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
      x = zeros(3)
      x[1:2] = mesh.vert_coords[:,vtx_arr[j],bndry_i.element]
      xi = x - ffd_box.origin
      for k = 1:3
        xi[k] = xi[k]/(ffd_box.box_bounds[2,k]-ffd_box.box_bounds[1,k])
      end
      err = xi - ffd_map.xi[itr][:,j,i]
      if norm(err) > 1e-14
        err_ctr += 1
      end
    end  # End for j = 1:length(vtx_arr)
  end    # End for i = 1:nfaces
end      # End for itr = 1:length(geomfaces)
if err_ctr > 0
  error("Counter > 0")
end
=#

#=
orig_cp_xyz = deepcopy(ffd_map.cp_xyz)
# Check if cp_xyz is updated correctly
err_ctr = 0
for i = 1:size(ffd_map.cp_xyz,2)
  for j = 1:size(ffd_map.cp_xyz,3)
    for k = 1:size(ffd_map.cp_xyz, 4)
      err = ffd_map.cp_xyz[1,i,j,k] - orig_cp_xyz[1,i,j,k]
      if norm(err - 0.02) > 1e-14
        err_ctr += 1
      end
    end
  end
end
if err_ctr > 0
  println("cp_xyz not updated")
end
=#

#=
err_ctr = 0
ctr2 = 0
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
      # fill!(x, 0.0)
      #println("original = $(orig_vert_coords[:,vtx_arr[j],bndry_i.element]), new = $(mesh.vert_coords[:,vtx_arr[j],bndry_i.element])")
      err = orig_vert_coords[:,vtx_arr[j],bndry_i.element] -
            mesh.vert_coords[:,vtx_arr[j],bndry_i.element] + [0.02,0.0]
      #println("err = $err")
      if norm(err) > 1e-14
        err_ctr += 1
      end
      ctr2 += 1
    end  # End for j = 1:length(vtx_arr)
  end    # End for i = 1:nfaces
end      # End for itr = 1:length(geomfaces)

println("ctr2 = $ctr2")
if err_ctr > 0
  println("Something wrong")
else
  println("err_ctr = $err_ctr")
end
=#

# Copy original mesh vertex coordinates
orig_vert_coords = deepcopy(mesh.vert_coords)

# Prep MeshWarping
volNodes = zeros(Tmsh, 3, mesh.numVert)
println("size mesh.vert_coords = $(size(mesh.vert_coords))")
for i = 1:mesh.numEl
  for j = 1:size(mesh.vert_coords,2)
    # Get the vertex numbering on the portion of mesh owned by the processor
    local_vertnum = mesh.element_vertnums[j,i]
    volNodes[1:2,local_vertnum] = mesh.vert_coords[:,j,i]# mesh.element_vertnums
    #println("volNodes[:, $local_vertnum] = $(volNodes[:,local_vertnum])")
  end
end

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
println("nwall_faces = $nwall_faces")
println("total surface vtx = ", sum(surfaceVtx))


wallCoords = zeros(3, sum(nwall_faces)*vtx_per_face)
# Populate wallCoords
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

# Warp parameters
nFaceConnectivity = sum(nwall_faces)*vtx_per_face
println("nFaceConnectivity = $nFaceConnectivity")
param = warpParam(nLocalNodes=mesh.numVert, nFaceConnLocal=nFaceConnectivity,
                  nFaceSizesLocal=sum(nwall_faces))
# MPI Information for warping
mpiVar = MPIInfo(warp_comm_world=comm_world, myID=my_rank, nProc=comm_size)

# Populate entries of param
param.aExp = 2.
param.bExp = 2.
param.LdefFact = 4.1
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
faceConn = collect(Int32, 0:nFaceConnectivity-1)
flatWarpSurfPts = reshape(wallCoords, 3*size(wallCoords,2))
faceSizes = mesh.dim*ones(Int32,sum(nwall_faces))
# warp the mesh

# Initialize warping
initializeWarping(param, mpiVar, symmetryPlanes, volNodes, surfaceVtx,
                      flatWarpSurfPts, faceConn, faceSizes)
# New Surface Coordinates

# Rotation matrix
theta = 10*pi/180  # Rotate wall coordinates by 10 degrees
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


# Translate control points along x & y by  5 units
# ffd_map.cp_xyz[1,:,:,:] += 0.02

evalSurface(ffd_map, mesh)
#=
# Check if evalSurface works correctly
err_ctr = 0
ctr2 = 0
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
      # fill!(x, 0.0)
      #println("original = $(orig_vert_coords[:,vtx_arr[j],bndry_i.element]), new = $(mesh.vert_coords[:,vtx_arr[j],bndry_i.element])")
      err = orig_vert_coords[:,vtx_arr[j],bndry_i.element] -
            mesh.vert_coords[:,vtx_arr[j],bndry_i.element] + [0.02,0.0]
      #println("err = $err")
      if norm(err) > 1e-14
        err_ctr += 1
      end
      ctr2 += 1
    end  # End for j = 1:length(vtx_arr)
  end    # End for i = 1:nfaces
end      # End for itr = 1:length(geomfaces)

println("ctr2 = $ctr2")
if err_ctr > 0
  println("Something wrong")
else
  println("All good. err_ctr = $err_ctr")
end
=#

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

# println("flatWarpSurfPts = \n",flatWarpSurfPts)

# Warp The mesh
warpMesh(param, volNodes, flatWarpSurfPts)

# Update all the mesh coordinates
println("type of volNodes = ", typeof(volNodes), " size = ", size(volNodes))
println("type of mesh.vert_coords = ", typeof(mesh.vert_coords), " size = ", size(mesh.vert_coords))
for i = 1:mesh.numEl
  for j = 1:size(mesh.vert_coords,2)
    # Get the vertex numbering on the portion of mesh owned by the processor
    local_vertnum = mesh.element_vertnums[j,i]
    mesh.vert_coords[:,j,i] = volNodes[1:2,local_vertnum] # mesh.element_vertnums
    #println("volNodes[:, $local_vertnum] = $(volNodes[:,local_vertnum])")
  end
end



for i = 1:mesh.numEl
  update_coords(mesh, i, mesh.vert_coords[:,:,i])
end
commit_coords(mesh)
PumiInterface.writeVtkFiles("rotation10", mesh.m_ptr)
