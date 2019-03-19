cmake_minimum_required(VERSION 3.14)

project(lib_mesh_simpl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

find_package(OpenMP)

add_library(${PROJECT_NAME}
            ecol.hpp
            ecol.cpp
            neighbor.hpp
            neighbor.cpp
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

if(OpenMP_Found)
    target_link_libraries(${PROJECT_NAME} openmp)
endif()

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Werror)
