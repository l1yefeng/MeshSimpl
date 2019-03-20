cmake_minimum_required(VERSION 3.7)
project(example LANGUAGES CXX)

set(OpenGL_GL_PREFERENCE GLVND)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(LIBIGL REQUIRED)

include_directories(../src)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME}
                      igl::core
                      igl::opengl_glfw
                      lib_mesh_simpl)
