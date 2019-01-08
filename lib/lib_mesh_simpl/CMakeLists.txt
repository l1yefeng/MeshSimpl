cmake_minimum_required(VERSION 3.13)
project(lib_mesh_simpl
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_library(lib_mesh_simpl SHARED
            library.h
            library.cpp
            util.h
            qem_heap.h)

# Force c++11
target_compile_features(lib_mesh_simpl PUBLIC cxx_std_11)
set_target_properties(lib_mesh_simpl PROPERTIES CXX_EXTENSIONS OFF)

# Related to tests. No panic if Catch2 cannot be found
find_package(Catch2 QUIET)
if(Catch2_FOUND)
    message(STATUS "Library Catch2 was found")
    add_executable(tests
                   test.cpp
                   util.h
                   library.h)
    target_link_libraries(tests Catch2::Catch2)
else()
    message(WARNING "Library Catch2 was not found -- Skipping building of tests")
endif()
