cmake_minimum_required(VERSION 3.21)
project(lab2_1 C)

SET(CMAKE_C_COMPILER gcc)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -fopenmp")

add_executable(lab2_1 lab2_1.c NCGM.c NCGM.h Matrix.c Matrix.h)
target_link_libraries(lab2_1 -lm)

#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})


#SET(CMAKE_C_COMPILER mpicc)
#SET(CMAKE_CXX_COMPILER mpicxx)
#target_link_libraries(lab2_1 -lm)
