# All 3D Pumi mesh tests go here

# MPI Declarations
if !MPI.Initialized()
  MPI.Init()
end

resize!(ARGS, 1)
ARGS[1] = "./input_vals_3d.jl"

opts = PDESolver.read_input(ARGS[1])
sbp, mesh, pmesh, Tsol, Tres, Tmsh, mesh_time = createMeshAndOperator(opts, 1)

facts("--- Checking FFD on 3D serial DG Pumi meshes ---") do

  context("Check control point manipulation with nonlinear mapping on the entire geometry") do
    
    ndim = mesh.dim
    order = [4,4,2]  # Order of B-splines in the 3 directions
    nControlPts = [4,4,2]
    offset = [0., 0., 0.5] # No offset in the X & Y direction

    # Create a mapping object using nonlinear mapping
    map = initializeFFD(mesh, sbp, order, nControlPts, offset, true)

    # Rigid body rotation
    theta = -20*pi/180  # Rotate wall coordinates by 10 degrees
    rotMat = [cos(theta) -sin(theta) 0
              sin(theta) cos(theta)  0
              0          0           1] # Rotation matrix
    # Rotate the control points
    for k = 1:map.nctl[3]
      for j = 1:map.nctl[2]
        for i = 1:map.nctl[1]
          map.cp_xyz[:,i,j,k] = rotMat*map.cp_xyz[:,i,j,k]
        end
      end
    end

    # Rigid body translation
    map.cp_xyz[1,:,:,:] += 0.2
    map.cp_xyz[2,:,:,:] += 0.3
    map.cp_xyz[3,:,:,:] += 0.5

    vertices = evalVolume(map, mesh)
    commitToPumi(map, mesh, sbp, vertices)

    fname = "./testvalues/translation_plus_rotation_DG_3D_tet8cube.dat"
    # f = open(fname, "w")
    # for i = 1:length(mesh.vert_coords)
    #   println(f, mesh.vert_coords[i])
    # end
    # close(f)

    test_values = readdlm(fname)
    for i = 1:length(test_values)
      err = abs(test_values[i] - mesh.vert_coords[i])
      @fact err --> less_than(1e-14)
    end

    
    writeVisFiles(mesh, "perturbed_3D_mesh")

  end # End context("Check FFD with nonlinear mapping on the entire geometry")


end # End facts("--- Checking FFD on 3D serial DG Pumi meshes ---")
