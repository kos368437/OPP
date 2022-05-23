cmake_minimum_required(VERSION 3.21)
project(lab2_2 C)

SET(CMAKE_C_COMPILER gcc)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lm -fopenmp")

add_executable(lab2_2 lab2_2.c NCGM.c NCGM.h Matrix.c Matrix.h)
target_link_libraries(lab2_2 -lm)
