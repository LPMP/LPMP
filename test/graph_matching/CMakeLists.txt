add_executable(koopmans_beckmann_test_input koopmans_beckmann_test_input.cpp)
target_link_libraries(koopmans_beckmann_test_input LPMP graph_matching_koopmans_beckmann_input)
add_test(koopmans_beckmann_test_input koopmans_beckmann_test_input) 

add_executable(test_graph_matching_input_export test_graph_matching_input_export.cpp)
target_link_libraries(test_graph_matching_input_export LPMP MRF_factors graph_matching_frank_wolfe)
add_test(test_graph_matching_input_export test_graph_matching_input_export)

add_test(NAME test_graph_matching_instance_python_construction
    COMMAND ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_graph_matching_instance_python_construction.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/graph_matching/
    )
set_tests_properties(test_graph_matching_instance_python_construction
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src/graph_matching:$ENV{PYTHONPATH}")

add_test(NAME test_graph_matching_instance_python_input_evaluation
    COMMAND ${PYTHON_EXECUTABLE}  ${CMAKE_CURRENT_SOURCE_DIR}/test_graph_matching_instance_python_input_evaluation.py
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/graph_matching/
    )
set_tests_properties(test_graph_matching_instance_python_input_evaluation
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src/graph_matching:$ENV{PYTHONPATH}")