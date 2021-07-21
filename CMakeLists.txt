cmake_minimum_required(VERSION 3.10)

project(karu)

set (CMAKE_CXX_STANDARD 11)

find_package(OpenCL REQUIRED)
find_package( Threads )
include_directories(src)
include_directories(${OpenCL_INCLUDE_DIRS})

set(SPARSE_KERNELS_PATH "src/algebra/sparse/lib/kernels/opencl")
set(COMPUTE_LIB_KERNELS_PATH "src/algebra/compute/lib/kernels/opencl")

macro(set_kernel_var filename variable)
    file(STRINGS ${filename} TMP)
    set(TMP2 "R\"(")
    foreach(LINE ${TMP})
        set(
            TMP2
            "${TMP2}
            ${LINE}"
        )
    endforeach()
    set(TMP2 "${TMP2})\"")
    set(${variable} ${TMP2})
endmacro()

set_kernel_var("${COMPUTE_LIB_KERNELS_PATH}/reduce_kernel.cl" REDUCE_CL_KERNEL_SRC)
set_kernel_var("${COMPUTE_LIB_KERNELS_PATH}/scan_kernel.cl" SCAN_CL_KERNEL_SRC)
set_kernel_var("${COMPUTE_LIB_KERNELS_PATH}/radix_sort_kernel.cl" RADIX_SORT_CL_KERNEL_SRC)
set_kernel_var("${SPARSE_KERNELS_PATH}/bsMV_kernel.cl" BSMAT_VEC_MULT_CL_KERNEL_SRC)
set_kernel_var("${SPARSE_KERNELS_PATH}/bsMM_kernel.cl" BSMAT_MAT_MULT_CL_KERNEL_SRC)

# Config LIb kernels
configure_file(src/algebra/compute/lib/kernels/kernels.template.in src/algebra/compute/lib/kernels/kernels.cpp)

# Config Matrix kernels
configure_file(src/algebra/sparse/lib/kernels/kernels.template.in src/algebra/sparse/lib/kernels/kernels.cpp)

add_library(
    karu
    src/algebra/compute/OpenCL.cpp
    src/algebra/compute/Context.cpp
    src/algebra/compute/Buffer.cpp
    src/algebra/compute/Kernel.cpp
    src/algebra/compute/Program.cpp
    src/algebra/compute/lib/Sort.cpp
    src/algebra/compute/lib/Scan.cpp
    src/algebra/compute/lib/Reduce.cpp
    src/algebra/compute/lib/Commom.cpp
    src/algebra/compute/lib/kernels/kernels.cpp
    src/algebra/vector/Vector.cpp
    src/algebra/matrix/MatrixData.cpp
    src/algebra/matrix/MatrixMultiplayer.cpp
    src/algebra/matrix/MatrixAdder.cpp
    src/algebra/matrix/MatrixSubtractor.cpp
    src/algebra/sparse/SparseMatrixData.cpp
    src/algebra/sparse/SparseMatrixMultiplayer.cpp
    src/algebra/sparse/lib/kernels/kernels.cpp
    src/algebra/sparse/lib/HashMap.cpp
)
    
target_include_directories(karu PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_link_libraries(karu OpenCL::OpenCL)
target_link_libraries(karu ${CMAKE_THREAD_LIBS_INIT})

target_compile_definitions(karu PRIVATE DEBUG=1)

enable_testing()
add_subdirectory(tests)
