cmake_minimum_required(VERSION 3.22)
project(Grobner)


set(CMAKE_CXX_STANDARD 14)
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

find_package(MPI REQUIRED)
include_directories("E:/mpi//sdk/Include")
link_directories("E:/mpi/sdk/Lib/x64")
link_libraries(msmpi)
add_executable(Grobner main.cpp Grobner_Matrix.h)
