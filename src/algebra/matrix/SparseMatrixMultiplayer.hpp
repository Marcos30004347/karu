#pragma once

#include "algebra/matrix/SparseMatrixData.hpp"

namespace karu {
namespace algebra {

class SparseMatrixMultiplayer {
	public:
	static void sparseMVMultiplyThreaded(SparseMatrixData* A, SparseMatrixData* x, SparseMatrixData* y);
	static void sparseMVMultiplyGPU(SparseMatrixData* A, SparseMatrixData* x, SparseMatrixData* y);

	static void sparseMMMultiplyGPU(SparseMatrixData* A, SparseMatrixData* B, SparseMatrixData* C);

	private:
	static void sparseMVMultiplyThreadHandler(void* data);
};

}
}
