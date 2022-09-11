cmake_minimum_required(VERSION 3.17)
project(Simplex)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(Simplex
        main.cpp
        simplex.cpp
        simplex.h)
