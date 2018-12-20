cmake_minimum_required(VERSION 3.13)
project(mesh_simpl)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(lib/lib_mesh_simpl)

add_executable(mesh_simpl
        main.cpp)

target_link_libraries(mesh_simpl
        lib_mesh_simpl)
