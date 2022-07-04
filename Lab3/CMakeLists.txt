cmake_minimum_required(VERSION 3.3)
project(lab_3 C)

add_executable(lab_3 lab3.c ParallelMultiplication.c ParallelMultiplication.h Matrix.c ../Matrix.h)

set(CMAKE_C_COMPILER /bin/mpicc)
target_link_libraries(lab_3 -lm)