cmake_minimum_required(VERSION 3.29)
project(newTransportationProblem)

set(CMAKE_CXX_STANDARD 20)

add_executable(newTransportationProblem main.cpp
        Matrix.cpp
        MethodOfVogel.cpp
        NorthWestCornerMethod.cpp)
