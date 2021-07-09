#pragma once

#include "algebra/matrix/BlockSparseMatrixData.hpp"

namespace karu {
namespace algebra {

class BlockSparseMatrixMultiplayer {
	public:
	static void sparseMVMultiplyThreaded(BlockSparseMatrixData* A, BlockSparseMatrixData* x, BlockSparseMatrixData* y);
	static void sparseMVMultiplyGPU(BlockSparseMatrixData* A, BlockSparseMatrixData* x, BlockSparseMatrixData* y);

	private:
	static void sparseMVMultiplyThreadHandler(void* data);
};

}
}
