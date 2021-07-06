cmake_minimum_required(VERSION 3.10)

include(CTest)

project(AlbebraComputeTests)
add_executable(AlbebraComputeTests algebra/core/compute.cpp)
target_link_libraries(AlbebraComputeTests karu)
add_test(NAME AlbebraComputeTests COMMAND AlbebraComputeTests)
