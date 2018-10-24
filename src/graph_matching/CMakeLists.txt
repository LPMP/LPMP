add_subdirectory(eval)

add_library(graph_matching_input graph_matching_input.cpp)
target_link_libraries(graph_matching_input LPMP)

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
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file})
   target_link_libraries( ${executable_file} graph_matching_input LPMP m stdc++ pthread )
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
    target_link_libraries( ${executable_file} LPMP graph_matching_input m stdc++ pthread FW-MAP arboricity)
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
add_executable(multigraph_matching multigraph_matching.cpp)
target_link_libraries(multigraph_matching multigraph_matching_input LPMP multicut_kernighan_lin transform_multigraph_matching)

add_executable(multigraph_matching_linear multigraph_matching_linear.cpp)
target_link_libraries(multigraph_matching_linear multigraph_matching_input LPMP)

add_executable(multigraph_matching_tightening multigraph_matching_tightening.cpp)
target_link_libraries(multigraph_matching_tightening multigraph_matching_input LPMP multicut_kernighan_lin transform_multigraph_matching)
