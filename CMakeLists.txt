cmake_minimum_required(VERSION 3.13)
project(mesh_simpl)

set(CMAKE_CXX_STANDARD 11)

# build lib_mesh_simpl beforehand which is itself a cmake project
add_subdirectory(lib/lib_mesh_simpl)

# find eigen library (required by libigl)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
# include igl headers
include_directories(lib/libigl/include/igl
                    lib/lib_mesh_simpl)

add_executable(mesh_simpl
               main.cpp
               measure.h
               write_obj.h)

target_compile_features(mesh_simpl PUBLIC cxx_std_11)
set_target_properties(mesh_simpl PROPERTIES CXX_EXTENSIONS OFF)

target_link_libraries(mesh_simpl
                      lib_mesh_simpl
                      Eigen3::Eigen)
