cmake_minimum_required(VERSION 3.10)

project(karu)

set (CMAKE_CXX_STANDARD 11)

find_package(OpenCL REQUIRED)

include_directories(src)
include_directories(${OpenCL_INCLUDE_DIRS})

add_library(
    karu
    src/algebra/core/compute/OpenCL.cpp
    src/algebra/core/compute/ComputeContext.cpp
    src/algebra/core/compute/ComputeBuffer.cpp
    src/algebra/core/compute/ComputeKernel.cpp
    src/algebra/core/compute/ComputeProgram.cpp
)
    
target_include_directories(karu PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_link_libraries(karu OpenCL::OpenCL)

enable_testing()
add_subdirectory(tests)
