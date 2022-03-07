find_package(OpenMP)

add_library(horizon_tracking_uai_input horizon_tracking_uai_input.cpp)
target_link_libraries(horizon_tracking_uai_input LPMP) 

add_executable(horizon_tracking_fwmap horizon_tracking_fwmap.cpp)
target_link_libraries(horizon_tracking_fwmap LPMP MRF_factors FW-MAP arboricity horizon_tracking_uai_input)
if(OPENMP_FOUND)
    target_link_libraries(horizon_tracking_fwmap OpenMP::OpenMP_CXX)
endif()

add_executable(horizon_tracking_subgradient horizon_tracking_subgradient.cpp)
target_link_libraries(horizon_tracking_subgradient LPMP MRF_factors arboricity horizon_tracking_uai_input)
if(OPENMP_FOUND)
    target_link_libraries(horizon_tracking_subgradient OpenMP::OpenMP_CXX)
endif()

add_executable(horizon_tracking_MST horizon_tracking_MST.cpp)
target_link_libraries(horizon_tracking_MST LPMP MRF_factors arboricity mrf_uai_input)

