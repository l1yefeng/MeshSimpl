cmake_minimum_required(VERSION 3.14)

project(lib_mesh_simpl LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# library
set(SRC src)

set(SOURCES
    ${SRC}/ecol.hpp
    ${SRC}/ecol.cpp
    ${SRC}/neighbor.hpp
    ${SRC}/neighbor.cpp
    ${SRC}/post_proc.hpp
    ${SRC}/post_proc.cpp
    ${SRC}/pre_proc.hpp
    ${SRC}/pre_proc.cpp
    ${SRC}/qem_heap.hpp
    ${SRC}/qem_heap.cpp
    ${SRC}/simplify.hpp
    ${SRC}/simplify.cpp
    ${SRC}/types.hpp
    ${SRC}/util.hpp
    ${SRC}/util.cpp)

add_library(${PROJECT_NAME} ${SOURCES})

target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Werror)

# example

add_subdirectory(example)
