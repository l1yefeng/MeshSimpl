cmake_minimum_required(VERSION 3.7)
project(MeshSimpl VERSION ${MESH_SIMPL_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_library(${PROJECT_NAME}
            edge.cpp
            edge.hpp
            erasable.hpp
            faces.cpp
            faces.hpp
            neighbor.hpp
            proc.cpp
            proc.hpp
            qemheap.cpp
            qemheap.hpp
            simplify.cpp
            simplify.hpp
            types.hpp
            util.cpp
            util.hpp nonmaniring.cpp nonmaniring.hpp)

target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_11)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall)
