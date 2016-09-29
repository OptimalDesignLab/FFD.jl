# Startup2.jl
# Merging mesh movement and FreeFormDeformation

push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, joinpath(Pkg.dir("PumiInterface"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("SummationByParts"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("MeshMovement"), "src"))
push!(LOAD_PATH, joinpath(Pkg.dir("PDESolver"), "src/Utils"))

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
# mesh_info = Int[sbp.numnodes, mesh.numEl, length(geom_faces)]

ffd_map = PumiMapping{Tmsh}(ndim, order, nControlPts, mesh, full_geom=false,
geom_faces=[5])
calcKnot(ffd_map)
# println("ffd_map.edge_knot = \n", ffd_map.edge_knot)

# Create Bounding box
offset = [0., 0., 0.5]
ffd_box = PumiBoundingBox{Tmsh}(ffd_map, mesh, sbp, offset)

# Control points
controlPoint(ffd_map, ffd_box)
