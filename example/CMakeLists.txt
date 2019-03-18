cmake_minimum_required(VERSION 3.14)

project(example LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  OFF)

find_package(LIBIGL REQUIRED QUIET)

include_directories(../src)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} igl::core lib_mesh_simpl)
