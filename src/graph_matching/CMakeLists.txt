add_subdirectory(eval)

add_library(graph_matching_instance  graph_matching_instance.cpp)
target_link_libraries(graph_matching_instance LPMP)

add_library(multigraph_matching_triplet_consistency_factor multigraph_matching_triplet_consistency_factor.cpp)
target_link_libraries(multigraph_matching_triplet_consistency_factor LPMP)

add_library(multigraph_matching_synchronization multigraph_matching_synchronization.cpp)
target_link_libraries(multigraph_matching_synchronization LPMP) 

add_library(multigraph_matching_factors INTERFACE)
target_link_libraries(multigraph_matching_factors INTERFACE MRF_factors)
target_link_libraries(multigraph_matching_factors INTERFACE multicut_instance)
target_link_libraries(multigraph_matching_factors INTERFACE correlation_clustering_instance)
target_link_libraries(multigraph_matching_factors INTERFACE multicut_kernighan_lin)
target_link_libraries(multigraph_matching_factors INTERFACE multigraph_matching_triplet_consistency_factor)
target_link_libraries(multigraph_matching_factors INTERFACE multigraph_matching_synchronization)
target_link_libraries(multigraph_matching_factors INTERFACE transform_multigraph_matching)

add_library(graph_matching_input graph_matching_input.cpp)
target_link_libraries(graph_matching_input LPMP)

add_library(graph_matching_koopmans_beckmann_input graph_matching_koopmans_beckmann_input.cpp)
target_link_libraries(graph_matching_koopmans_beckmann_input LPMP)

add_library(multigraph_matching_input multigraph_matching_input.cpp)
target_link_libraries(multigraph_matching_input LPMP graph_matching_input)

SET(SOURCE_FILES
   graph_matching_mp_left.cpp graph_matching_mp_right.cpp graph_matching_mp_both_sides.cpp 
   graph_matching_mcf_left.cpp graph_matching_mcf_right.cpp graph_matching_mcf_both_sides.cpp 
   graph_matching_gm_left.cpp graph_matching_gm_right.cpp 
   #graph_matching_mp_uai.cpp graph_matching_mcf_uai.cpp graph_matching_gm_uai.cpp
   graph_matching_hungarian_bp_left.cpp graph_matching_hungarian_bp_right.cpp graph_matching_hungarian_bp_both_sides.cpp 
   #hungarian_bp_uai.cpp

   # with frustrated cycle tightening
   graph_matching_mp_left_tightening.cpp graph_matching_mp_right_tightening.cpp graph_matching_mp_both_sides_tightening.cpp 
   graph_matching_mcf_left_tightening.cpp graph_matching_mcf_right_tightening.cpp graph_matching_mcf_both_sides_tightening.cpp 
   graph_matching_gm_left_tightening.cpp graph_matching_gm_right_tightening.cpp 
   #graph_matching_mp_uai_tightening.cpp graph_matching_mcf_uai_tightening.cpp graph_matching_gm_uai_tightening.cpp
   graph_matching_hungarian_bp_left_tightening.cpp graph_matching_hungarian_bp_right_tightening.cpp graph_matching_hungarian_bp_both_sides_tightening.cpp 
   #hungarian_bp_uai_tightening.cpp

   graph_matching_mp_both_sides_inter_quadratic_message graph_matching_mp_both_sides_inter_quadratic_message_tightening
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file})
   target_link_libraries( ${executable_file} MRF_factors graph_matching_input LPMP)
endforeach( source_file ${SOURCE_FILES} )

# with proximal bundle
SET(PROXIMAL_BUNDLE_SOLVERS
  graph_matching_mcf_proximal_bundle_left.cpp
  graph_matching_gm_proximal_bundle_left.cpp 
  #graph_matching_local_subproblems_proximal_bundle_left.cpp
  )

foreach( source_file ${PROXIMAL_BUNDLE_SOLVERS} )
    string( REPLACE ".cpp" "" executable_file ${source_file} )
    add_executable( ${executable_file} ${source_file} )
    target_link_libraries( ${executable_file} LPMP MRF_factors graph_matching_input m stdc++ pthread FW-MAP arboricity)
endforeach( source_file ${SOURCE_FILES} )

# with conic bundle
#add_executable(graph_matching_mcf_conic_bundle graph_matching_mcf_conic_bundle.cpp ${headers} ${sources})
#target_link_libraries( graph_matching_mcf_conic_bundle m stdc++ pthread lgl SVM arboricity cb)

if(WITH_GUROBI)
   SET(GUROBI_FILES
      # lp solver
      graph_matching_mp_left_tightening_gurobi.cpp graph_matching_mp_right_tightening_gurobi.cpp graph_matching_mp_both_sides_tightening_gurobi.cpp 
      graph_matching_mcf_left_tightening_gurobi.cpp graph_matching_mcf_right_tightening_gurobi.cpp graph_matching_mcf_both_sides_tightening_gurobi.cpp 
      graph_matching_gm_left_tightening_gurobi.cpp graph_matching_gm_right_tightening_gurobi.cpp 
      #graph_matching_mp_uai_tightening_gurobi.cpp graph_matching_mcf_uai_tightening_gurobi.cpp graph_matching_gm_uai_tightening_gurobi.cpp
      )

   foreach( source_file ${GUROBI_FILES} )
      string( REPLACE ".cpp" "" executable_file ${source_file} )
      add_executable( ${executable_file} ${source_file} ${headers} ${sources})
      target_link_libraries( ${executable_file} m stdc++ pthread ${GUROBI_LIBRARIES})

      target_compile_definitions(${executable_file} PUBLIC -DUSE_GUROBI)
      target_include_directories(${executable_file} PUBLIC ${GUROBI_INCLUDE_DIR})
   endforeach( source_file ${GUROBI_FILES} )
endif()

# multigraph matching
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
   target_link_libraries( ${executable_file} multigraph_matching_input LPMP multigraph_matching_factors)
endforeach( source_file ${MGM_SOURCE_FILES} )

add_executable(graph_matching_fw graph_matching_fw.cpp)
target_link_libraries(graph_matching_fw LPMP graph_matching_input)

add_executable(graph_matching_frank_wolfe_sparse_text_input graph_matching_frank_wolfe_sparse_text_input.cpp)
target_link_libraries(graph_matching_frank_wolfe_sparse_text_input LPMP graph_matching_input)

add_executable(graph_matching_frank_wolfe_koopmans_beckmann_text_input graph_matching_frank_wolfe_koopmans_beckmann_text_input.cpp)
target_link_libraries(graph_matching_frank_wolfe_koopmans_beckmann_text_input LPMP graph_matching_koopmans_beckmann_input)
