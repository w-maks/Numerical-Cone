cmake_minimum_required(VERSION 3.29)
project(Numerical_Cone)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(Numerical_Cone
    rk4.cpp
    rk4.h
        cone.cpp
        cone.h
        main.cpp
)
