cmake_minimum_required(VERSION 3.7)
project(lib_mesh_simpl LANGUAGES CXX)

option(LIB_MESH_SIMPL_EXAMPLE "Build example executable" ON)
option(LIB_MESH_SIMPL_TEST    "Build unit tests"         ON)

add_subdirectory(src)

# example
if(LIB_MESH_SIMPL_EXAMPLE)
    add_subdirectory(example)
endif()

# test
if(LIB_MESH_SIMPL_TEST)
    add_subdirectory(tests)
endif()
