cmake_minimum_required(VERSION 3.24)

project(roots)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
add_definitions("-Wall -O2")

enable_testing()

add_subdirectory(src)
add_subdirectory(test)