add_executable(dimacs_max_flow_input_test dimacs_max_flow_input_test.cpp)
target_link_libraries(dimacs_max_flow_input_test LPMP dimacs_max_flow_input)
add_test(dimacs_max_flow_input_test dimacs_max_flow_input_test)

add_executable(generate_export_read_binary_MRF_instance_test generate_export_read_binary_MRF_instance_test.cpp)
target_link_libraries(generate_export_read_binary_MRF_instance_test LPMP mrf_uai_input)
add_test(generate_export_read_binary_MRF_instance_test generate_export_read_binary_MRF_instance_test)

add_executable(binary_pairwise_potential_test binary_pairwise_potential_test.cpp)
target_link_libraries(binary_pairwise_potential_test LPMP)
add_test(binary_pairwise_potential_test binary_pairwise_potential_test)

add_executable(transform_max_flow_test transform_max_flow_test.cpp)
target_link_libraries(transform_max_flow_test LPMP)
add_test(transform_max_flow_test transform_max_flow_test)

add_executable(potts_factor potts_factor.cpp)
target_link_libraries(potts_factor LPMP m stdc++ pthread lingeling)
add_test(potts_factor potts_factor)

add_executable(simplex simplex.cpp)
target_link_libraries(simplex LPMP m stdc++ pthread lingeling)
add_test(simplex simplex)

add_executable(simplex_marginalization simplex_marginalization.cpp)
target_link_libraries(simplex_marginalization LPMP m stdc++ pthread lingeling)
add_test(simplex_marginalization simplex_marginalization)

add_executable(cycle_inequalities cycle_inequalities.cpp)
target_link_libraries(cycle_inequalities LPMP m stdc++ pthread lingeling)
add_test(cycle_inequalities cycle_inequalities)

add_executable(mrf_tree_decomposition mrf_tree_decomposition.cpp)
target_link_libraries(mrf_tree_decomposition LPMP m stdc++ pthread lingeling FW-MAP arboricity cb)
add_test(mrf_tree_decomposition mrf_tree_decomposition)

add_executable(mrf_tree mrf_tree.cpp)
target_link_libraries(mrf_tree LPMP m stdc++ pthread lingeling arboricity)
add_test(mrf_tree mrf_tree)

add_executable(mrf_read_uai_file mrf_read_uai_file.cpp)
target_link_libraries(mrf_read_uai_file LPMP m stdc++ pthread lingeling mrf_uai_input)
add_test(mrf_read_uai_file mrf_read_uai_file)

add_executable(mrf_sat_rounding mrf_sat_rounding.cpp)
target_link_libraries(mrf_sat_rounding LPMP m stdc++ pthread lingeling)
add_test(mrf_sat_rounding mrf_sat_rounding)

add_executable(mrf_tightening mrf_tightening.cpp)
target_link_libraries(mrf_tightening LPMP m stdc++ pthread lingeling)
add_test(mrf_tightening mrf_tightening)

if(WITH_GUROBI)
  add_executable(combiLP combiLP.cpp)
  target_link_libraries(combiLP LPMP m stdc++ pthread GUROBI HDF5 lingeling) #somehow currently needs to be linked against. Remove this!
  add_test(combiLP combiLP)

  add_executable(dd_ilp_gurobi dd_ilp_gurobi.cpp)
  target_link_libraries(dd_ilp_gurobi LPMP m stdc++ pthread GUROBI HDF5 lingeling) 
  add_test(dd_ilp_gurobi dd_ilp_gurobi)
endif()
