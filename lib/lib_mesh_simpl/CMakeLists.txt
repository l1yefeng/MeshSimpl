cmake_minimum_required(VERSION 3.13)
project(lib_mesh_simpl)

set(CMAKE_CXX_STANDARD 14)

add_library(lib_mesh_simpl SHARED
        library.cpp
        library.h)
