#pragma once

// #include "algebra/vector/Vector.hpp"
#include "algebra/sparse/SparseMatrixData.hpp"
#include "algebra/matrix/MatrixData.hpp"

namespace karu::algebra {

class SparseMatrixMultiplayer {
public:
	static void sparseMVMultiplyCPU(SparseMatrixData* A, MatrixData* x, MatrixData* y);
	static void sparseMVMultiplyGPU(SparseMatrixData* A, MatrixData* x, MatrixData* y);

	// static void sparseMMMultiplyGPU(const SparseMatrixData* const A, const SparseMatrixData* const B, SparseMatrixData* C);

private:
	static void sparseMVMultiplyThreadHandler(void* data);
};

}
