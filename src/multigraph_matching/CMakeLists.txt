add_subdirectory(eval)

add_library(multigraph_matching_triplet_consistency_factor multigraph_matching_triplet_consistency_factor.cpp)
target_link_libraries(multigraph_matching_triplet_consistency_factor LPMP)

add_library(multigraph_matching_synchronization multigraph_matching_synchronization.cpp)
target_link_libraries(multigraph_matching_synchronization LPMP) 

add_library(multigraph_matching_instance multigraph_matching_instance.cpp)
target_link_libraries(multigraph_matching_instance LPMP graph_matching_instance)

add_library(multigraph_matching_factors INTERFACE)
target_link_libraries(multigraph_matching_factors INTERFACE MRF_factors)
target_link_libraries(multigraph_matching_factors INTERFACE multicut_instance)
target_link_libraries(multigraph_matching_factors INTERFACE correlation_clustering_instance)
target_link_libraries(multigraph_matching_factors INTERFACE multicut_kernighan_lin)
target_link_libraries(multigraph_matching_factors INTERFACE multigraph_matching_triplet_consistency_factor)
target_link_libraries(multigraph_matching_factors INTERFACE multigraph_matching_synchronization)
target_link_libraries(multigraph_matching_factors INTERFACE transform_multigraph_matching)

add_library(multigraph_matching_input multigraph_matching_input.cpp)
target_link_libraries(multigraph_matching_input LPMP graph_matching_input)

add_executable(multigraph_matching_universe_fw multigraph_matching_universe_fw.cpp)
target_link_libraries(multigraph_matching_universe_fw LPMP multigraph_matching_input) 

SET(MGM_SOURCE_FILES
      multigraph_matching_mp.cpp
      multigraph_matching_mcf.cpp
      multigraph_matching_mcf_inter_quadratic_messages.cpp
      multigraph_matching_tightening_mp.cpp
      multigraph_matching_tightening_mcf.cpp
      multigraph_matching_tightening_mcf_inter_quadratic_message.cpp
      evaluate_multigraph_matching.cpp
      )

foreach( source_file ${MGM_SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file})
   target_link_libraries( ${executable_file} multigraph_matching_input LPMP multigraph_matching_factors graph_matching_frank_wolfe)
endforeach( source_file ${MGM_SOURCE_FILES} )

pybind11_add_module(multigraph_matching_py multigraph_matching_python_binding.cpp)
target_link_libraries(multigraph_matching_py PRIVATE LPMP graph_matching_py_helper graph_matching_frank_wolfe multigraph_matching_factors multigraph_matching_input MRF_factors)
