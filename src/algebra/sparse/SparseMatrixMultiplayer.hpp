#pragma once

#include "algebra/vector/Vector.hpp"
#include "algebra/sparse/SparseMatrixData.hpp"

namespace karu::algebra {

class SparseMatrixMultiplayer {
	public:
	static void sparseMVMultiplyCPU(SparseMatrixData* A, Vector* x, Vector* y);
	static void sparseMVMultiplyGPU(SparseMatrixData* A, Vector* x, Vector* y);

	static void sparseMMMultiplyGPU(SparseMatrixData* A, SparseMatrixData* B, SparseMatrixData* C);

	private:
	static void sparseMVMultiplyThreadHandler(void* data);
};

}
