add_subdirectory(mrf)
add_subdirectory(graph_matching)
add_subdirectory(max_cut)
add_subdirectory(horizon_tracking)
add_subdirectory(multicut)
add_subdirectory(bdd)

add_executable(test_message_passing_schedule test_message_passing_schedule.cpp)
target_link_libraries(test_message_passing_schedule LPMP m stdc++)
add_test(test_message_passing_schedule test_message_passing_schedule) 

add_executable(test_weight_computation test_weight_computation.cpp)
target_link_libraries(test_weight_computation LPMP m stdc++)
add_test(test_weight_computation test_weight_computation) 

add_executable(test_topological_sort test_topological_sort.cpp)
target_link_libraries(test_topological_sort LPMP m stdc++)
add_test(test_topological_sort test_topological_sort) 

add_executable(test_two_dimensional_variable_array test_two_dimensional_variable_array.cpp)
target_link_libraries(test_two_dimensional_variable_array  LPMP m stdc++)
add_test(test_two_dimensional_variable_array  test_two_dimensional_variable_array) 

add_executable(test_three_dimensional_variable_array test_three_dimensional_variable_array.cpp)
target_link_libraries(test_three_dimensional_variable_array  LPMP m stdc++)
add_test(test_three_dimensional_variable_array  test_three_dimensional_variable_array) 

add_executable(test_reparametrization_recognition test_reparametrization_recognition.cpp)
target_link_libraries(test_reparametrization_recognition  LPMP m stdc++)
add_test(test_reparametrization_recognition  test_reparametrization_recognition) 

add_executable(vector vector.cpp)
target_link_libraries( vector LPMP m stdc++ pthread )
add_test( vector vector )

add_executable(serialization serialization.cpp)
target_link_libraries( serialization LPMP m stdc++ pthread )
add_test( serialization serialization )

add_executable(test_model test_model.cpp)
target_link_libraries(test_model LPMP DD_ILP lingeling)
add_test( test_model test_model )

add_executable(test_FWMAP test_FWMAP.cpp)
target_link_libraries(test_FWMAP LPMP FW-MAP lingeling)
add_test(test_FWMAP test_FWMAP)

add_executable(test_conic_bundle test_conic_bundle.cpp)
target_link_libraries(test_conic_bundle cb LPMP lingeling)
add_test(test_conic_bundle test_conic_bundle)

add_executable(graph_test graph_test.cpp)
target_link_libraries(graph_test LPMP)
add_test(graph_test graph_test)

add_executable(test_commutation_matrix test_commutation_matrix.cpp)
target_link_libraries(test_commutation_matrix LPMP)
add_test(test_commutation_matrix test_commutation_matrix)

add_executable(test_polynomial_functions test_polynomial_functions.cpp)
target_link_libraries(test_polynomial_functions LPMP)
add_test(test_polynomial_function test_polynomial_functions)

add_executable(test_cuthill_mckee test_cuthill_mckee.cpp)
target_link_libraries(test_cuthill_mckee LPMP)
add_test(test_cuthill_mckee test_cuthill_mckee)

add_executable(test_bfs_ordering test_bfs_ordering.cpp)
target_link_libraries(test_bfs_ordering LPMP)
add_test(test_bfs_ordering test_bfs_ordering)
