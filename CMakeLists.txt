project(LPMP)

cmake_minimum_required(VERSION 2.8.12)

# C++11
add_compile_options(-std=c++17)

# compiler options
add_definitions(-march=native)

add_library(LPMP INTERFACE)
target_include_directories(LPMP INTERFACE external/qpbo/include)
target_include_directories(LPMP INTERFACE external/maxflow/include)
target_include_directories(LPMP INTERFACE include)
target_include_directories(LPMP INTERFACE external/PEGTL)
target_include_directories(LPMP INTERFACE external/TCLAP/include)
target_include_directories(LPMP INTERFACE external/libsimdpp)
target_include_directories(LPMP INTERFACE external/meta/include)
target_include_directories(LPMP INTERFACE external/BCFW-Bundle/include)
target_include_directories(LPMP INTERFACE external/ConicBundle/include)
target_include_directories(LPMP INTERFACE external/arboricity)
target_include_directories(LPMP INTERFACE external/MCF-SSP)
target_include_directories(LPMP INTERFACE external/andres-graph/include)

# for opengm #
target_include_directories(LPMP INTERFACE external/opengm/include)

find_package(OpenMP)
target_link_libraries(LPMP INTERFACE OpenMP::OpenMP_CXX)

find_package(HDF5 REQUIRED CXX)
add_library(HDF5 INTERFACE)
target_include_directories(HDF5 INTERFACE ${HDF5_INCLUDE_DIR})
target_link_libraries(HDF5 INTERFACE ${HDF5_LIBRARIES})
target_link_libraries(HDF5 INTERFACE ${HDF5_CXX_LIBRARIES})
target_compile_definitions(HDF5 INTERFACE ${HDF5_DEFINITIONS})
message(STATUS ${HDF5_LIBRARIES})
message(STATUS ${HDF5_INCLUDE_DIR})
##############

enable_testing()

add_subdirectory(external/qpbo)
add_subdirectory(external/BCFW-Bundle) 
add_subdirectory(external/DD_ILP)
add_subdirectory(external/ConicBundle)
add_subdirectory(external/arboricity)

target_link_libraries(LPMP INTERFACE DD_ILP)

add_subdirectory(src)
enable_testing()
add_subdirectory(test)

option(BUILD_DOC "Build documentation" OFF)
if(BUILD_DOC)
    find_package(Doxygen REQUIRED)
configure_file(Doxyfile.in Doxyfile @ONLY)
    add_custom_target(docs-lpmp
            COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM)
    if(NOT TARGET docs)
    add_custom_target(docs ALL)
    endif()
    add_dependencies(docs docs-lpmp)
endif()
