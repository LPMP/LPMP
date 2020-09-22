
add_library(ldp_instance ldp_instance.cxx)
target_link_libraries(ldp_instance LPMP  pybind11::module)
 
#add_library(lifted_disjoint_paths_input lifted_disjoint_paths_input.cpp)
#target_link_libraries(lifted_disjoint_paths_input LPMP ldp_instance pybind11::module)

add_executable(lifted_disjoint_paths_text_input  lifted_disjoint_paths_text_input.cpp)
target_link_libraries(lifted_disjoint_paths_text_input ldp_instance LPMP pybind11::module)

pybind11_add_module(ldpPython ldp_python.cxx)
target_link_libraries(ldpPython PRIVATE ldp_instance LPMP pybind11::module)   


