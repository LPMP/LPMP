add_library(ldp_directed_graph ldp_directed_graph.cxx)
target_link_libraries(ldp_directed_graph LPMP pybind11::module)


add_library(ldp_instance ldp_instance.cxx )

target_link_libraries(ldp_instance LPMP  pybind11::module ldp_directed_graph)

add_library(ldp_two_layer_graph ldp_two_layer_graph.cxx)
target_link_libraries(ldp_two_layer_graph LPMP)

add_library(ldp_cut_factor ldp_cut_factor.cxx)
target_link_libraries(ldp_cut_factor LPMP ldp_two_layer_graph)

add_library(ldp_path_factor ldp_path_factor.cxx)
target_link_libraries(ldp_path_factor LPMP pybind11::module)


add_library(ldp_batch_process ldp_batch_process.cxx)
target_link_libraries(ldp_batch_process LPMP pybind11::module)
 

#add_library(lifted_disjoint_paths_input lifted_disjoint_paths_input.cpp)
#target_link_libraries(lifted_disjoint_paths_input LPMP ldp_instance pybind11::module)

add_executable(lifted_disjoint_paths_text_input  lifted_disjoint_paths_text_input.cpp)

target_link_libraries(lifted_disjoint_paths_text_input ldp_instance ldp_cut_factor ldp_path_factor ldp_directed_graph ldp_batch_process LPMP pybind11::module)

pybind11_add_module(ldpMessagePassingPy ldp_python.cxx)

target_link_libraries(ldpMessagePassingPy PRIVATE ldp_instance  ldp_cut_factor ldp_path_factor ldp_directed_graph  ldp_batch_process LPMP pybind11::module)

configure_file(solveFromFiles.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromFiles.py)
#configure_file(solveFromFilesForBatch.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromFilesForBatch.py)
configure_file(solveFromVectors.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromVectors.py)
configure_file(solveFromVectorsTwoGraphs.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromVectorsTwoGraphs.py)
configure_file(solveInBatches.py ${CMAKE_CURRENT_BINARY_DIR}/solveInBatches.py)
#configure_file(solveFromExamples.py ${CMAKE_CURRENT_BINARY_DIR}/solveFromExamples.py)



