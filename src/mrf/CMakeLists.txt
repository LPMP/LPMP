add_library(MRF_factors pairwise_simplex_factor.cpp pairwise_Potts_factor.cpp ternary_simplex_factor)
target_link_libraries(MRF_factors LPMP) 

add_library(dimacs_max_flow_input dimacs_max_flow_input.cpp)
target_link_libraries(dimacs_max_flow_input LPMP)

add_executable(convert_opengm_mrf_to_uai convert_opengm_mrf_to_uai.cpp)
target_link_libraries(convert_opengm_mrf_to_uai LPMP HDF5)

add_executable(convert_dimacs_qpbo_max_flow_to_mrf_uai convert_dimacs_qpbo_max_flow_to_mrf_uai.cpp)
target_link_libraries(convert_dimacs_qpbo_max_flow_to_mrf_uai LPMP dimacs_max_flow_input)

add_executable(convert_dimacs_graph_cut_max_flow_to_mrf_uai convert_dimacs_graph_cut_max_flow_to_mrf_uai.cpp)
target_link_libraries(convert_dimacs_graph_cut_max_flow_to_mrf_uai LPMP dimacs_max_flow_input) 

add_library(mrf_uai_input mrf_uai_input.cpp)
target_link_libraries(mrf_uai_input LPMP) 

add_library(mrf_opengm_input mrf_opengm_input.cpp)
target_link_libraries(mrf_opengm_input LPMP HDF5) 

SET(UAI_SOURCE_FILES
   mplp_uai.cpp
   srmp_uai.cpp
   srmp_uai_tightening.cpp
   qpbo_uai.cpp
   )

foreach( source_file ${UAI_SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} )
   target_link_libraries( ${executable_file} LPMP MRF_factors mrf_uai_input ) 
endforeach( source_file ${SOURCE_FILES} )

SET(OPENGM_SOURCE_FILES
   mplp_opengm.cpp
   srmp_opengm.cpp
   srmp_opengm_tightening.cpp
   )

foreach( source_file ${OPENGM_SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
   add_executable( ${executable_file} ${source_file} ${headers} ${sources})
   target_link_libraries( ${executable_file} LPMP MRF_factors mrf_opengm_input ) 
endforeach( source_file ${SOURCE_FILES} )

add_executable(FWMAP_uai FWMAP_uai.cpp)
target_link_libraries(FWMAP_uai LPMP MRF_factors FW-MAP arboricity mrf_uai_input)

add_executable(FWMAP_opengm FWMAP_opengm.cpp)
target_link_libraries(FWMAP_opengm HDF5 LPMP MRF_factors FW-MAP arboricity mrf_opengm_input) 

add_executable(conic_bundle_uai conic_bundle_uai.cpp)
target_link_libraries(conic_bundle_uai LPMP MRF_factors cb arboricity mrf_uai_input)

add_executable(conic_bundle_opengm conic_bundle_opengm.cpp)
target_link_libraries(conic_bundle_opengm HDF5 LPMP MRF_factors cb arboricity mrf_opengm_input)

add_executable(qpbo_dimacs qpbo_dimacs.cpp)
target_link_libraries(qpbo_dimacs LPMP dimacs_max_flow_input)
