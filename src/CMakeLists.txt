cmake_minimum_required(VERSION 3.7)
project(lib_mesh_simpl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)

add_library(${PROJECT_NAME}
            ecol.hpp
            ecol.cpp
            edge.hpp
            neighbor.hpp
            post_proc.hpp
            post_proc.cpp
            pre_proc.hpp
            pre_proc.cpp
            qem_heap.hpp
            qem_heap.cpp
            simplify.hpp
            simplify.cpp
            types.hpp
            util.hpp
            util.cpp)

if(OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic)
