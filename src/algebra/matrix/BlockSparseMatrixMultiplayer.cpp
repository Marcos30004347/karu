#include <assert.h>
#include <iostream>
#include <vector>
#include <thread>

#include "algebra/core/compute/Compute.hpp"
#include "algebra/matrix/kernels/kernels.hpp"
#include "algebra/matrix/BlockSparseMatrixMultiplayer.hpp"

using namespace karu;
using namespace algebra;

struct SparseMVMultiplyData
{
	u64 idx;

	BlockSparseMatrixData* A;
	BlockSparseMatrixData* x;
	BlockSparseMatrixData* y;

	bool A_transposed;
};

void BlockSparseMatrixMultiplayer::sparseMVMultiplyThreadHandler(void* data)
{
	SparseMVMultiplyData* args = (SparseMVMultiplyData*)data;

	u64 target_block_row = args->idx;

	BlockSparseMatrixData* A = args->A;
	BlockSparseMatrixData* x = args->x;
	BlockSparseMatrixData* y = args->y;

	u64 first_block = A->bcsr_row_ptr[target_block_row];
	u64 last_block = A->bcsr_row_ptr[target_block_row + 1];

	u64 bs = A->bcsr_block_width * A->bcsr_block_heigth;
	u64 v_bs = x->bcsr_block_width * x->bcsr_block_heigth;

	std::vector<f32> local_out = std::vector<f32>(bs, 0);

	for(u64 block = first_block; block < last_block; block++)
	{

		u64 target_block_col = A->bcsr_col_idx[block];

		for(u64 c = 0; c<A->bcsr_block_width; c++)
		{
			u64 vec_this_col = x->bcsr_data[target_block_col + c];
			for(u64 r = 0; r<A->bcsr_block_heigth; r++)
			{
				local_out[r] += A->bcsr_data[block*bs + c*A->bcsr_block_heigth + r] * vec_this_col;
			}
		}
	}

	for(u64 r=0; r<v_bs; r++)
	{
		y->bcsr_data[target_block_row*v_bs + r] += local_out[r];
	}
}

void BlockSparseMatrixMultiplayer::sparseMVMultiplyThreaded(BlockSparseMatrixData* A, BlockSparseMatrixData* x, BlockSparseMatrixData* y)
{
	u64 n = A->bcsr_lines/A->bcsr_block_heigth;

	std::thread t[n];
	SparseMVMultiplyData data[n];

	for(u64 i=0; i<n; i++)
	{
		data[i].idx = i;
		data[i].A = A;
		data[i].x = x;
		data[i].y = y;

		t[i] = std::thread(BlockSparseMatrixMultiplayer::sparseMVMultiplyThreadHandler, &data[i]);
	}

	for(u64 i=0; i<n; i++)
		t[i].join();
}

void BlockSparseMatrixMultiplayer::sparseMVMultiplyGPU(BlockSparseMatrixData* A, BlockSparseMatrixData* x, BlockSparseMatrixData* y)
{
	compute::Buffer col_idx = compute::Buffer(A->bcsr_col_idx.data(), sizeof(u64)*A->bcsr_col_idx.size(), compute::Buffer::READ_ONLY, false);
	compute::Buffer row_ptr = compute::Buffer(A->bcsr_row_ptr.data(), sizeof(u64)*A->bcsr_row_ptr.size(), compute::Buffer::READ_ONLY, false);
	compute::Buffer A_buffe	= compute::Buffer(A->bcsr_data.data(), 		sizeof(f32)*A->bcsr_data.size(), 		compute::Buffer::READ_ONLY, false);
	compute::Buffer x_buffe = compute::Buffer(x->bcsr_data.data(), 		sizeof(f32)*x->bcsr_data.size(),	 	compute::Buffer::READ_ONLY, false);

	compute::Buffer y_buffe = compute::Buffer(sizeof(f32)*y->bcsr_data.size(), compute::Buffer::WRITE_ONLY, compute::Buffer::MEM_GPU);

	const unsigned int block_dim = A->bcsr_block_heigth * A->bcsr_block_width;

	bsMV_kernel->setKernelArgument(0, sizeof(unsigned long int), &A->bcsr_block_heigth);
	bsMV_kernel->setKernelArgument(1, sizeof(unsigned long int), &A->bcsr_block_width);
	bsMV_kernel->setKernelArgument(2, sizeof(unsigned long int), &x->bcsr_block_heigth);
	bsMV_kernel->setKernelArgument(3, sizeof(unsigned long int), &x->bcsr_block_width);
	bsMV_kernel->setKernelArgument(4, BUFFER_ARG_SIZE, col_idx.computeUnitRef());
	bsMV_kernel->setKernelArgument(5, BUFFER_ARG_SIZE, row_ptr.computeUnitRef());
	bsMV_kernel->setKernelArgument(6, BUFFER_ARG_SIZE, A_buffe.computeUnitRef());
	bsMV_kernel->setKernelArgument(7, BUFFER_ARG_SIZE, x_buffe.computeUnitRef());
	bsMV_kernel->setKernelArgument(8, BUFFER_ARG_SIZE, y_buffe.computeUnitRef());

	// shared memory
	bsMV_kernel->setKernelArgument(9, block_dim * sizeof(float), nullptr);

	bsMV_kernel->enqueue({ A->bcsr_data.size() }, { block_dim });

	y_buffe.toLogicUnit();

	f32* data = (f32*)y_buffe.logicUnitRef();

	for(u64 i=0; i<y->bcsr_lines; i++)
	{
		y->bcsr_data[i] = data[i];
		// std::cout << data[i] << " ";
	}

	// std::cout << std::endl;
	// std::cout << std::endl;
}