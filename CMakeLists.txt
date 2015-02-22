cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(graph_contraction)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}) #"${CMAKE_SOURCE_DIR}/../MLUtils/cmake/"

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Build directory" FORCE)

#find_package(Eigen REQUIRED)
#include_directories(${EIGEN_INCLUDE_DIRS})
#find_package(Boost COMPONENTS filesystem system REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost COMPONENTS python REQUIRED)

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
include_directories("include/")

add_executable(example src/example.cpp)
#add_library(graphcontraction src/graphcontraction.cpp)
PYTHON_ADD_MODULE(graphcontraction src/graphcontraction.cpp)
target_link_libraries(graphcontraction ${Boost_LIBRARIES})
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-Wall -g -Werror -std=c++11 ${CMAKE_CXX_FLAGS}")
