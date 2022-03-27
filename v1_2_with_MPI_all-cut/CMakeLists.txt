cmake_minimum_required(VERSION 3.3)
project(lab1_2 C)

add_executable(lab1_2 lab1_2.c NCGMP.c NCGMP.h ../Matrix.c ../Matrix.h)

set(CMAKE_C_COMPILER /usr/bin/mpicc)