cmake_minimum_required(VERSION 3.10)

include(CTest)

project(AlgebraComputeTests)
add_executable(AlgebraComputeTests algebra/compute/Compute.cpp)
target_link_libraries(AlgebraComputeTests karu)
add_test(NAME AlgebraComputeTests COMMAND AlgebraComputeTests)

project(AlgebraComputeScanTests)
add_executable(AlgebraComputeScanTests algebra/compute/Scan.cpp)
target_link_libraries(AlgebraComputeScanTests karu)
add_test(NAME AlgebraComputeScanTests COMMAND AlgebraComputeScanTests)

project(AlgebraComputeReduceTests)
add_executable(AlgebraComputeReduceTests algebra/compute/Reduce.cpp)
target_link_libraries(AlgebraComputeReduceTests karu)
add_test(NAME AlgebraComputeReduceTests COMMAND AlgebraComputeReduceTests)

project(AlgebraComputeSortTests)
add_executable(AlgebraComputeSortTests algebra/compute/Sort.cpp)
target_link_libraries(AlgebraComputeSortTests karu)
add_test(NAME AlgebraComputeSortTests COMMAND AlgebraComputeSortTests)

project(AlgebraMatrixDataTests)
add_executable(AlgebraMatrixDataTests algebra/matrix/MatrixData.cpp)
target_link_libraries(AlgebraMatrixDataTests karu)
add_test(NAME AlgebraMatrixDataTests COMMAND AlgebraMatrixDataTests)

project(AlgebraTensorTests)
add_executable(AlgebraTensorTests algebra/tensor/Tensor.cpp)
target_link_libraries(AlgebraTensorTests karu)
add_test(NAME AlgebraTensorTests COMMAND AlgebraTensorTests)

project(AlgebraMatrixMultiplayerTests)
add_executable(AlgebraMatrixMultiplayerTests algebra/matrix/MatrixMultiplayer.cpp)
target_link_libraries(AlgebraMatrixMultiplayerTests karu)
add_test(NAME AlgebraMatrixMultiplayerTests COMMAND AlgebraMatrixMultiplayerTests)

project(AlgebraMatrixTests)
add_executable(AlgebraMatrixTests algebra/matrix/Matrix.cpp)
target_link_libraries(AlgebraMatrixTests karu)
add_test(NAME AlgebraMatrixTests COMMAND AlgebraMatrixTests)

project(AlgebraSparseMatrixTests)
add_executable(AlgebraSparseMatrixTests algebra/sparse/SpMatrix.cpp)
target_link_libraries(AlgebraSparseMatrixTests karu)
add_test(NAME AlgebraSparseMatrixTests COMMAND AlgebraSparseMatrixTests)

project(AlgebraSparseMatrixDataTests)
add_executable(AlgebraSparseMatrixDataTests algebra/sparse/SparseMatrixData.cpp)
target_link_libraries(AlgebraSparseMatrixDataTests karu)
add_test(NAME AlgebraSparseMatrixDataTests COMMAND AlgebraSparseMatrixDataTests)

project(AlgebraSparseMatrixMultiplayerTests)
add_executable(AlgebraSparseMatrixMultiplayerTests algebra/sparse/SparseMatrixMultiplayer.cpp)
target_link_libraries(AlgebraSparseMatrixMultiplayerTests karu)
add_test(NAME AlgebraSparseMatrixMultiplayerTests COMMAND AlgebraSparseMatrixMultiplayerTests)

project(LinearRotationTests)
add_executable(LinearRotationTests algebra/linear/Rotation.cpp)
target_link_libraries(LinearRotationTests karu)
add_test(NAME LinearRotationTests COMMAND LinearRotationTests)


project(LinearConjugateGradientTests)
add_executable(LinearConjugateGradientTests algebra/linear/ConjugateGradient.cpp)
target_link_libraries(LinearConjugateGradientTests karu)
add_test(NAME LinearConjugateGradientTests COMMAND LinearConjugateGradientTests)

project(BundleAdjustmentTests)
add_executable(BundleAdjustmentTests bundle/BundleAdjustment.cpp)
target_link_libraries(BundleAdjustmentTests karu)
add_test(NAME BundleAdjustmentTests COMMAND BundleAdjustmentTests)

project(PolynomialTests)
add_executable(PolynomialTests algebra/polynomial/Polynomial.cpp)
target_link_libraries(PolynomialTests karu)
add_test(NAME PolynomialTests COMMAND PolynomialTests)

project(HouseholderTests)
add_executable(HouseholderTests algebra/SVD/Householder.cpp)
target_link_libraries(HouseholderTests karu)
add_test(NAME HouseholderTests COMMAND HouseholderTests)

project(GolubKahanTests)
add_executable(GolubKahanTests algebra/SVD/GolubKahan.cpp)
target_link_libraries(GolubKahanTests karu)
add_test(NAME GolubKahanTests COMMAND GolubKahanTests)

project(SVDTests)
add_executable(SVDTests algebra/SVD/SVD.cpp)
target_link_libraries(SVDTests karu)
add_test(NAME SVDTests COMMAND SVDTests)

project(QRTests)
add_executable(QRTests algebra/SVD/QR.cpp)
target_link_libraries(QRTests karu)
add_test(NAME QRTests COMMAND QRTests)

project(GaussianTests)
add_executable(GaussianTests gaussian/test1.cpp)
target_link_libraries(GaussianTests karu)
add_test(NAME GaussianTests COMMAND GaussianTests)

project(CameraTests)
add_executable(CameraTests camera/Camera.cpp)
target_link_libraries(CameraTests karu)
add_test(NAME CameraTests COMMAND CameraTests)
