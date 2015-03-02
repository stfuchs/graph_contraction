cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(graph_contraction)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Build directory" FORCE)
set(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Build directory" FORCE)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

include_directories("include/")

add_executable(example src/example.cpp)
add_executable(qh_example src/qh_example.cpp)

#add_library(eigenmap src/eigenmap.cpp)
#add_library(libeigenmap SHARED py/libeigenmap.cpp)
#set_target_properties(libeigenmap PROPERTIES PREFIX "")
#target_link_libraries(libeigenmap eigenmap)
add_library(graphcontraction SHARED py/graphcontraction.cpp)
set_target_properties(graphcontraction PROPERTIES PREFIX "")

#PYTHON_ADD_MODULE(graphcontraction src/graphcontraction.cpp)
#target_link_libraries(graphcontraction ${Boost_LIBRARIES})
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wall -g -std=c++11 ${CMAKE_CXX_FLAGS}")
