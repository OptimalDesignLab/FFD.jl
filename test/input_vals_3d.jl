# this user supplied file creates a dictionary of arguments
# if a key is repeated, the last use of the key is used
# it is a little bit dangerous letting the user run arbitrary code
# as part of the solver
# now that this file is read inside a function, it is better encapsulated

arg_dict = Dict{ASCIIString, Any}(
"run_type" => 5,
"jac_method" => 1,
"jac_type" => 2,
"order" => 1,
"use_DG" => true,
"dimensions" => 3,
"Flux_name" => "RoeFlux",
"IC_name" => "ICExp",
"variable_type" => :conservative,
"numBC" => 1,
"BC1" => [ 0, 1, 2, 3],
"BC1_name" => "ExpBC",
#"BC2" => [4, 10],
#"BC2_name" => "noPenetrationBC",
"delta_t" => 0.005,
"t_max" => 500.000,
"smb_name" => "SRCMESHES/tet8cube.smb",
"dmg_name" => ".null",
"res_abstol" => 1e-12,
"res_reltol" => 1e-10,
"step_tol" => 1e-10,
"itermax" => 10000,
"tau_type" => 3,
"edgestab_gamma" => -0.9,
"use_filter" => false,
"use_dissipation" => false,
"dissipation_name" => "damp1",
"dissipation_const" => 12.00,
#"calc_error" => true,
#"calc_error_infname" => "solution_final.dat",
#"writeq" => true,
#"perturb_ic" => true,
#"perturb_mag" => 0.001,
#"write_sparsity" => true,
#"write_jac" => true,
#"write_edge_vertnums" => true,
#"write_face_vertnums" => true,
#"write_qic" => true,
#"writeboundary" => true,
#"write_res" => true,
"write_finalsolution" => true,
"write_finalresidual" => true,
"write_counts" => true,
"write_rhs" => false,
"do_postproc" => true,
"exact_soln_func" => "ICIsentropicVortex",
"solve" => false,
)
