cmake_minimum_required(VERSION 3.7)
project(MeshSimpl VERSION ${MESH_SIMPL_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_library(${PROJECT_NAME}
        connectivity.hpp
        ecol.cpp
        ecol.hpp
        edge.cpp
        edge.hpp
        marker.hpp
        neighbor.hpp
        post_proc.cpp
        post_proc.hpp
        pre_proc.cpp
        pre_proc.hpp
        qem_heap.cpp
        qem_heap.hpp
        simplify.cpp
        simplify.hpp
        types.hpp
        util.cpp
        util.hpp)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)

target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
