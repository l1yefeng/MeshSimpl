cmake_minimum_required(VERSION 3.13)
project(lib_mesh_simpl
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

add_library(${PROJECT_NAME}
            ecol.h
            ecol.cpp
            util.h
            util.cpp
            qem_heap.h
            qem_heap.cpp
            simplify.h
            simplify.cpp
            pre_proc.h
            pre_proc.cpp
            post_proc.h
            post_proc.cpp
            neighbor.cpp
            neighbor.h)

# Force c++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -pedantic -Werror)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_EXTENSIONS OFF)

# Related to tests. No panic if Catch2 cannot be found
find_package(Catch2 QUIET)
if(Catch2_FOUND)
    message(STATUS "Library Catch2 was found")
    add_executable(tests test.cpp ${SOURCES})
    target_link_libraries(tests Catch2::Catch2)
else()
    message(WARNING "Library Catch2 was not found -- Skipping building of tests")
endif()
