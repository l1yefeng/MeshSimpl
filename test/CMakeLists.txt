cmake_minimum_required(VERSION 3.14)

project(test LANGUAGES CXX)

include_directories(../src)

find_package(Catch2)
if(Catch2_FOUND)
    add_executable(${PROJECT_NAME} test.cpp)
    target_link_libraries(${PROJECT_NAME} lib_mesh_simpl Catch2::Catch2)
endif()