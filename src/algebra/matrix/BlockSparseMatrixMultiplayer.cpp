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

void BlockSparseMatrixMultiplayer::sparseMMMultiplyGPU(BlockSparseMatrixData* A, BlockSparseMatrixData* B, BlockSparseMatrixData* C)
{

	// Get overlaping blocks that will be multiplyed
	u64 bs = A->bcsr_block_heigth;
	for(u32 i=0; i<A->lines()/A->bcsr_block_heigth; i++)
	{
		u32 A_first_block_row = A->bcsr_row_ptr[i];
		u32 A_last_block_row = A->bcsr_row_ptr[i+1];

		for(u32 A_block = A_first_block_row; A_block < A_last_block_row; A_block++)
		{
			
			for(u32 j=0; j<B->lines()/B->bcsr_block_heigth; j++)
			{
				u32 B_first_block = B->bcsr_row_ptr[j];
				u32 B_last_block = B->bcsr_row_ptr[j+1];
				
				for(u32 B_block=B_first_block; B_block<B_last_block; B_block++)
				{
					i64 B_col = B->bcsr_col_idx[B_block];
					i64 A_col = A->bcsr_col_idx[A_block];
				
					u32 cond0 = i*(u32)A->bcsr_block_heigth - j*(u32)B->bcsr_block_heigth;
					u32 cond1 = A_col - j*(u32)B->bcsr_block_heigth;
					u32 dist0 = cond0 > 0 ? cond0 : -cond0;
					u32 dist1 = cond1 > 0 ? cond1 : -cond1;

					if(dist1 < A->bcsr_block_width)
					{
						// A_block * B_block
						std::cout << A_block << " * " << B_block << std::endl;
					}
				}
			}
		}
	}
	// for(u64 line=0; line<block_lines; line++)
	// {
	// 	std::cout << "line: " << line << std::endl;

	// 	u64 first_block = A->bcsr_row_ptr[line];
	// 	u64 last_block = A->bcsr_row_ptr[line+1];

	// 	for(u64 block = first_block; block<last_block; block++)
	// 	{
	// 		u64 A_column = A->bcsr_col_idx[block];
	// 		u64 A_block  = A_column/bs;
	
	// 		u64 B_block  = B->bcsr_row_ptr[A_block]/bs;
			
	// 		std::cout << "T: " << A_column/bs << std::endl;
	// 		std::cout << "X: " << B->bcsr_row_ptr[A_column/bs] << std::endl;
	// 		std::cout << "W: " << B->bcsr_col_idx[B->bcsr_row_ptr[A_column/bs]] << std::endl;
	// 		// std::cout << "A_block: " << A_block << std::endl;

	// 		// if(B_block == A_block)
	// 		// {
	// 		// 	std::cout << "Ab(" << A_block << ") * Bb(" << B_block  << ")" << std::endl;
	// 		// }
	// 	}
	// 	std::cout << std::endl;
	// }

	// compute::Buffer A_col_idx = compute::Buffer(A->bcsr_col_idx.data(), sizeof(u64)*A->bcsr_col_idx.size(), compute::Buffer::READ_ONLY, false);
	// compute::Buffer A_row_ptr = compute::Buffer(A->bcsr_row_ptr.data(), sizeof(u64)*A->bcsr_row_ptr.size(), compute::Buffer::READ_ONLY, false);
	// compute::Buffer B_col_idx = compute::Buffer(B->bcsr_col_idx.data(), sizeof(u64)*B->bcsr_col_idx.size(), compute::Buffer::READ_ONLY, false);
	// compute::Buffer B_row_ptr = compute::Buffer(B->bcsr_row_ptr.data(), sizeof(u64)*B->bcsr_row_ptr.size(), compute::Buffer::READ_ONLY, false);

	// compute::Buffer A_buffe	= compute::Buffer(A->bcsr_data.data(), 		sizeof(f32)*A->bcsr_data.size(), compute::Buffer::READ_WRITE, false);
	// compute::Buffer B_buffe = compute::Buffer(B->bcsr_data.data(), 		sizeof(f32)*B->bcsr_data.size(), compute::Buffer::READ_WRITE, false);
	// compute::Buffer C_buffe = compute::Buffer(sizeof(f32)*C->bcsr_data.size(), compute::Buffer::WRITE_ONLY, compute::Buffer::MEM_GPU);

	// const unsigned int block_dim = A->bcsr_block_heigth * A->bcsr_block_width;

	// bsMM_kernel->setKernelArgument(0, sizeof(unsigned long int), &A->bcsr_lines);
	// bsMM_kernel->setKernelArgument(1, sizeof(unsigned long int), &A->bcsr_columns);
	// bsMM_kernel->setKernelArgument(2, sizeof(unsigned long int), &B->bcsr_lines);
	// bsMM_kernel->setKernelArgument(3, sizeof(unsigned long int), &B->bcsr_columns);
	// bsMM_kernel->setKernelArgument(4, sizeof(unsigned long int), &A->bcsr_block_heigth);
	// bsMM_kernel->setKernelArgument(5, sizeof(unsigned long int), &A->bcsr_block_width);
	// bsMM_kernel->setKernelArgument(6, BUFFER_ARG_SIZE, A_col_idx.computeUnitRef());
	// bsMM_kernel->setKernelArgument(7, BUFFER_ARG_SIZE, A_row_ptr.computeUnitRef());
	// bsMM_kernel->setKernelArgument(8, BUFFER_ARG_SIZE, B_col_idx.computeUnitRef());
	// bsMM_kernel->setKernelArgument(9, BUFFER_ARG_SIZE, B_row_ptr.computeUnitRef());
	// bsMM_kernel->setKernelArgument(10, BUFFER_ARG_SIZE, A_buffe.computeUnitRef());
	// bsMM_kernel->setKernelArgument(11, BUFFER_ARG_SIZE, B_buffe.computeUnitRef());
	// bsMM_kernel->setKernelArgument(12, BUFFER_ARG_SIZE, C_buffe.computeUnitRef());
	// bsMM_kernel->setKernelArgument(13, block_dim * sizeof(float), nullptr);

	// bsMM_kernel->enqueue({ A->bcsr_data.size() }, { block_dim });

	// // C_buffe.toLogicUnit();
	// A_buffe.toLogicUnit();
	// B_buffe.toLogicUnit();

	// f32* A_data = (f32*)A_buffe.logicUnitRef();
	// f32* B_data = (f32*)B_buffe.logicUnitRef();
	// // f32* data = (f32*)C_buffe.logicUnitRef();

	// // A->bcsr_data.clear();
	// // B->bcsr_data.clear();

	// for(int i=0; i<A->bcsr_data.size(); i++)
	// 		A->bcsr_data[i] = A_data[i];

	// for(int i=0; i<B->bcsr_data.size(); i++)
	// 		B->bcsr_data[i] = B_data[i];

	// for(int i=0; i<A->lines(); i++)
	// {
	// 	for(int j=0; j<A->columns(); j++)
	// 	{
	// 		std::cout << A->get(i, j) << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	// std::cout << std::endl;
	// for(int i=0; i<B->lines(); i++)
	// {
	// 	for(int j=0; j<B->columns(); j++)
	// 	{
	// 		std::cout << B->get(i, j) << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

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
}
