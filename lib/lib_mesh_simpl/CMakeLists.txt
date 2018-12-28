cmake_minimum_required(VERSION 3.13)
project(lib_mesh_simpl)

set(CMAKE_CXX_STANDARD 14)

if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
elseif (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
endif ()

add_library(lib_mesh_simpl SHARED
            library.h
            library.cpp
            util.h)
