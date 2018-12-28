cmake_minimum_required(VERSION 3.13)
project(mesh_simpl)

set(CMAKE_CXX_STANDARD 14)

# build lib_mesh_simpl beforehand which is itself a cmake project
add_subdirectory(lib/lib_mesh_simpl)

# find libraries
find_package(Threads)

# find eigen library (required by libigl)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
# include igl headers
include_directories(lib/libigl/include/igl
                    lib/lib_mesh_simpl)

add_executable(mesh_simpl
               main.cpp
               measure.h)

target_link_libraries(mesh_simpl
                      lib_mesh_simpl
                      Eigen3::Eigen
                      ${CMAKE_THREAD_LIBS_INIT})
