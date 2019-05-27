cmake_minimum_required(VERSION 3.7)
project(mesh_simpl LANGUAGES CXX)

set(MESH_SIMPL_VERSION 1.0.0)

option(LIB_MESH_SIMPL_EXAMPLE "Build example executable" ON)
option(LIB_MESH_SIMPL_BENCH "Build example executable" OFF)

add_subdirectory(src)


# example
if(LIB_MESH_SIMPL_EXAMPLE)
  set(TARGET_EXAMPLE "example")
  add_executable(${TARGET_EXAMPLE} example/main.cpp)
  target_include_directories(${TARGET_EXAMPLE} PRIVATE src)
  target_link_libraries(${TARGET_EXAMPLE} MeshSimpl)
endif()

# bench
if(LIB_MESH_SIMPL_BENCH)
  add_subdirectory(bench)
endif()
