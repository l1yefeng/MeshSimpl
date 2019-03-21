cmake_minimum_required(VERSION 3.7)
project(lib_mesh_simpl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

set(LIB_MESH_SIMPL_OPENMP "Enable OpenMP parallelization" OFF)

add_library(${PROJECT_NAME}
            connectivity.hpp
            ecol.cpp
            ecol.hpp
            edge.hpp
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

if(LIB_MESH_SIMPL_OPENMP)
    find_package(OpenMP)
    if(OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic")
