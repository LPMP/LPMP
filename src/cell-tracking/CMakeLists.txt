add_library(cell_tracking_input cell_tracking_input.cpp)
target_link_libraries(cell_tracking_input LPMP)

add_executable(convert_cell_tracking_to_lp convert_cell_tracking_to_lp.cpp)
target_link_libraries(convert_cell_tracking_to_lp cell_tracking_input LPMP)

SET(SOURCE_FILES
  cell_tracking_mother_machine.cpp 
  cell_tracking_with_division_distance.cpp
  cell_tracking.cpp 
  cell_tracking_flow.cpp 
  cell_tracking_division_distance_duplicate_edges.cpp 
  cell_tracking_duplicate_edges.cpp 
  cell_tracking_fine_decomposition.cpp
  )

foreach( source_file ${SOURCE_FILES} )
   string( REPLACE ".cpp" "" executable_file ${source_file} )
#   add_executable( ${executable_file} ${source_file})
#   target_link_libraries( ${executable_file} LP_MP lingeling m stdc++ pthread cell_tracking_input )
endforeach( source_file ${SOURCE_FILES} )

#add_executable( comparison_plots  comparison_plots.cpp)
#target_link_libraries( comparison_plots LP_MP lingeling m stdc++ pthread sqlite3 )

