cmake_minimum_required(VERSION 3.10)

project(karu)

set (CMAKE_CXX_STANDARD 11)

find_package(OpenCL REQUIRED)
find_package( Threads )
include_directories(src)
include_directories(${OpenCL_INCLUDE_DIRS})

set(MATRIX_KERNELS_PATH "src/algebra/matrix/kernels/opencl")

macro(set_kernel_var filename variable)
    file(STRINGS ${filename} TMP)
    set(TMP2 " ")
    foreach(LINE ${TMP})
        set(
            TMP2
            "${TMP2}
            \"${LINE}\\n\""
        )
    endforeach()
    set(${variable} ${TMP2})
endmacro()

set_kernel_var("${MATRIX_KERNELS_PATH}/bsMV_kernel.cl" BSMAT_VEC_MULT_CL_KERNEL_SRC)

configure_file(src/algebra/matrix/kernels/kernels.template.in src/algebra/matrix/kernels/kernels.cpp)

add_library(
    karu
    src/algebra/core/compute/OpenCL.cpp
    src/algebra/core/compute/ComputeContext.cpp
    src/algebra/core/compute/ComputeBuffer.cpp
    src/algebra/core/compute/ComputeKernel.cpp
    src/algebra/core/compute/ComputeProgram.cpp
    src/algebra/matrix/MatrixData.cpp
    src/algebra/matrix/MatrixMultiplayer.cpp
    src/algebra/matrix/MatrixAdder.cpp
    src/algebra/matrix/MatrixSubtractor.cpp
    src/algebra/matrix/BlockSparseMatrixData.cpp
    src/algebra/matrix/BlockSparseMatrixMultiplayer.cpp
    src/algebra/matrix/kernels/kernels.cpp
)
    
target_include_directories(karu PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_link_libraries(karu OpenCL::OpenCL)
target_link_libraries(karu ${CMAKE_THREAD_LIBS_INIT})

target_compile_definitions(karu PRIVATE DEBUG=1)

enable_testing()
add_subdirectory(tests)
