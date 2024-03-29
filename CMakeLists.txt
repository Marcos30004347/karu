cmake_minimum_required(VERSION 3.10)

project(karu)

set (CMAKE_CXX_STANDARD 11)


add_subdirectory(vendor/glfw)
find_package(OpenGL REQUIRED)
find_package(OpenCV REQUIRED)

# include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIRS})

# link_libraries(${GLEW_LIBRARIES})
link_libraries(${OPENGL_LIBRARIES})
link_libraries(${OPENGL_gl_LIBRARY})
link_libraries(${OPENGL_glu_LIBRARY})


set(GLAD_DIR "vendor/glad")
add_library("glad" "${GLAD_DIR}/src/glad.c")
target_include_directories("glad" PRIVATE "${GLAD_DIR}/include")

find_package(OpenCL REQUIRED)
find_package( Threads )
include_directories(src)
include_directories(${OpenCL_INCLUDE_DIRS})

set(SPARSE_KERNELS_PATH "src/algebra/sparse/lib/kernels/opencl")
set(COMPUTE_LIB_KERNELS_PATH "src/algebra/compute/lib/kernels/opencl")
set(SIFT_KERNELS_PATH "src/sift/kernels/opencl")


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

set_kernel_var("${SIFT_KERNELS_PATH}/gaussian_blur.cl" BLUR_CL_KERNEL_SRC)
set_kernel_var("${SIFT_KERNELS_PATH}/resize.cl" RESIZE_KERNEL_SRC)

# Config LIb kernels
configure_file(src/algebra/compute/lib/kernels/kernels.template.in src/algebra/compute/lib/kernels/kernels.cpp)

# Config Matrix kernels
configure_file(src/algebra/sparse/lib/kernels/kernels.template.in src/algebra/sparse/lib/kernels/kernels.cpp)

# Config SIFT kernels
configure_file(src/sift/kernels/kernels.template.in src/sift/kernels/kernels.cpp)

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
    src/camera/codegen/CameraModelBundle.cpp
    src/camera/codegen/CameraModel.cpp
    src/bundle/codegen/Homography.cpp
    src/bundle/codegen/Estimation.cpp
    src/algebra/vector/Vector.cpp
    src/algebra/polynomial/Polynomial.cpp
    src/algebra/matrix/Matrix.cpp
    src/algebra/matrix/MatrixEchelonForm.cpp
    src/algebra/matrix/MatrixNullSpace.cpp
    src/algebra/matrix/MatrixData.cpp
    src/algebra/matrix/MatrixMultiplayer.cpp
    src/algebra/matrix/MatrixAdder.cpp
    src/algebra/matrix/MatrixSubtractor.cpp
    src/algebra/matrix/MatrixTransposer.cpp
    src/algebra/matrix/MatrixDivider.cpp
    src/algebra/matrix/MatrixLU.cpp
    src/algebra/sparse/SpMatrix.cpp
    src/algebra/sparse/SparseMatrixData.cpp
    src/algebra/sparse/SparseMatrixMultiplayer.cpp
    src/algebra/sparse/lib/kernels/kernels.cpp
    src/algebra/sparse/lib/HashMap.cpp
    src/algebra/linear/Linear.cpp
    src/renderer/Renderer.cpp
    
    src/sift/GaussianBlur.cpp
    src/sift/Resize.cpp
    src/sift/kernels/kernels.cpp
    src/sift/Sift.cpp
    src/sift/Matcher.cpp
    
    src/point-cloud/PointCloud.cpp

    src/image/Image.cpp
)

# set_property(TARGET karu PROPERTY
#     CXX_STANDARD 11
#     CXX_STANDARD_REQUIRED ON
# )

target_include_directories(karu PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_include_directories(karu PUBLIC "${GLAD_DIR}/include")

target_link_libraries(karu OpenCL::OpenCL)
target_link_libraries(karu glfw)
target_link_libraries(karu ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(karu glad)
target_link_libraries(karu ${OpenCV_LIBS} )

target_compile_definitions(karu PRIVATE DEBUG=1)

add_subdirectory(samples)

enable_testing()
add_subdirectory(tests)
