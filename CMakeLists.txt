cmake_minimum_required(VERSION 3.10)

project(karu)

set (CMAKE_CXX_STANDARD 11)

add_library(
    karu
    karu/compute.cpp
    karu/matrix.cpp
)

target_include_directories(karu PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

