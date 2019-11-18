add_library(graph_matching_instance  graph_matching_instance.cpp)
target_link_libraries(graph_matching_instance LPMP)

add_library(graph_matching_input graph_matching_input.cpp)
target_link_libraries(graph_matching_input LPMP)

add_library(graph_matching_koopmans_beckmann_input graph_matching_koopmans_beckmann_input.cpp)
target_link_libraries(graph_matching_koopmans_beckmann_input LPMP)

add_library(graph_matching_frank_wolfe graph_matching_frank_wolfe.cpp)
target_link_libraries(graph_matching_frank_wolfe LPMP)

#add_library(graph_matching_FMC_mp graph_matching_FMC_mp.cpp)
#target_link_libraries(graph_matching_FMC_mp LPMP MRF_factors)

#add_executable(lpmp_graph_matching lpmp_graph_matching.cpp)
#target_link_libraries(lpmp_graph_matching LPMP graph_matching_input MRF_factors)

SET(SOURCE_FILES
   graph_matching_mp.cpp 
   graph_matching_gm.cpp 
   graph_matching_hungarian_bp.cpp 

   # with frustrated cycle tightening
   graph_matching_mp_tightening.cpp
   graph_matching_gm_tightening.cpp
   graph_matching_hungarian_bp_tightening.cpp 

   # with inter-quadratic messages
   graph_matching_mp_inter_quadratic_message 
   graph_matching_mp_inter_quadratic_message_tightening
   )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file})
   target_link_libraries( ${executable_file} MRF_factors graph_matching_input graph_matching_frank_wolfe LPMP)
endforeach( source_file ${SOURCE_FILES} )

# with proximal bundle
SET(PROXIMAL_BUNDLE_SOLVERS
  graph_matching_mcf_proximal_bundle.cpp
  graph_matching_gm_proximal_bundle.cpp 
  )

foreach( source_file ${PROXIMAL_BUNDLE_SOLVERS} )
    string( REPLACE ".cpp" "" executable_file ${source_file} )
    add_executable( ${executable_file} ${source_file} )
    target_link_libraries( ${executable_file} LPMP MRF_factors graph_matching_input graph_matching_frank_wolfe FW-MAP arboricity)
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

add_executable(graph_matching_frank_wolfe_text_input graph_matching_frank_wolfe_text_input.cpp)
target_link_libraries(graph_matching_frank_wolfe_text_input LPMP graph_matching_input graph_matching_frank_wolfe) 

add_executable(graph_matching_frank_wolfe_koopmans_beckmann_text_input graph_matching_frank_wolfe_koopmans_beckmann_text_input.cpp)
target_link_libraries(graph_matching_frank_wolfe_koopmans_beckmann_text_input LPMP graph_matching_koopmans_beckmann_input)

pybind11_add_module(graph_matching_python_binding graph_matching_python_binding.cpp)
target_link_libraries(graph_matching_python_binding PRIVATE LPMP graph_matching_input MRF_factors graph_matching_frank_wolfe)
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
CONFIGURE_FILE(setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
