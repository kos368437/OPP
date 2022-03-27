cmake_minimum_required(VERSION 3.3)
project(lab1_1 C)

add_executable(lab1_1 lab1_1.c NCGMP.c NCGMP.h ../Matrix.c ../Matrix.h)

set(CMAKE_C_COMPILER /usr/bin/mpicc)
