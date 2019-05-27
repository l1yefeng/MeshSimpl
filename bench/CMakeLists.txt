
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)

add_subdirectory(benchmark)

add_executable(bench
               bench.cpp)

target_include_directories(bench
                           PRIVATE
                           ../src)

target_link_libraries(bench
                      PRIVATE
                      benchmark
                      MeshSimpl)