cmake_minimum_required(VERSION 2.8)
project(2dv CXX)

set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall")

find_package(OpenMP)
if (OPENMP_FOUND)
  message("OPENMP FOUND")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_subdirectory(restructure)
add_executable(
  ${PROJECT_NAME}
  2dv.cpp
  _class_and_variables.cpp
  _class.cpp
  _parameters.cpp
  output.cpp
  force.cpp
  ODE_solver.cpp
)
target_link_libraries(${PROJECT_NAME} m restructure)
