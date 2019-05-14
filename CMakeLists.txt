cmake_minimum_required(VERSION 3.7)
project(mesh_simpl LANGUAGES CXX)

set(MESH_SIMPL_VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)

option(LIB_MESH_SIMPL_EXAMPLE "Build example executable" ON)
option(LIB_MESH_SIMPL_TEST "Build unit tests" OFF)

add_subdirectory(src)


# example
if(LIB_MESH_SIMPL_EXAMPLE)
  set(TARGET_EXAMPLE "example")
  add_executable(${TARGET_EXAMPLE} example/main.cpp)
  target_include_directories(${TARGET_EXAMPLE} PRIVATE src)
  target_link_libraries(${TARGET_EXAMPLE} MeshSimpl)
endif()

# test
if(LIB_MESH_SIMPL_TEST)
  find_package(Catch2)
  if(Catch2_FOUND)
    set(TARGET_TEST "tests")
    add_executable(${TARGET_TEST} tests/test.cpp)
    target_include_directories(${TARGET_TEST} PRIVATE src)
    target_link_libraries(${TARGET_TEST} MeshSimpl Catch2::Catch2)
  endif()
endif()
