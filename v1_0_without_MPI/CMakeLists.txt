cmake_minimum_required(VERSION 3.21)
project(lab1_0 C)

add_executable(lab1_0 lab1_0.c NCGM.c NCGM.h ../Matrix.c ../Matrix.h)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

