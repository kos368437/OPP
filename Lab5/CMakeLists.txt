cmake_minimum_required(VERSION 3.22)
project(Lab_5 C)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)

set(CMAKE_C_STANDARD 11)

add_executable(Lab_5 main.c )
target_link_libraries(Lab_5 -lpthread -lm)
