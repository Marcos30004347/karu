#include <assert.h>
#include <iostream>
#include <vector>
#include <thread>

#include "algebra/sparse/lib/kernels/kernels.hpp"

#include "algebra/sparse/lib/Sort.hpp"
#include "algebra/sparse/lib/HashMap.hpp"
#include "algebra/compute/Compute.hpp"
#include "algebra/sparse/SparseMatrixMultiplayer.hpp"

using namespace karu;
using namespace algebra;

struct SparseMVMultiplyData
{
	u64 idx;

	SparseMatrixData* A;
	Vector* x;
	Vector* y;

	bool A_transposed;
};

void SparseMatrixMultiplayer::sparseMVMultiplyThreadHandler(void* data)
{
	SparseMVMultiplyData* args = (SparseMVMultiplyData*)data;

	u64 target_block_row = args->idx;

	SparseMatrixData* A = args->A;

	Vector* x = args->x;
	Vector* y = args->y;

	f32* x_data = x->data();
	f32* y_data = y->data();

	u64 first_block = A->bcsr_row_ptr[target_block_row];
	u64 last_block = A->bcsr_row_ptr[target_block_row + 1];

	u64 bs = A->bcsr_block_width * A->bcsr_block_heigth;

	std::vector<f32> local_out = std::vector<f32>(bs, 0);
	for(u64 block = first_block; block < last_block; block++)
	{
		u64 target_block_col = A->bcsr_col_idx[block];

		for(u64 c = 0; c<A->bcsr_block_width; c++)
		{
			f32 vec_this_col = x_data[target_block_col + c];
			for(u64 r = 0; r<A->bcsr_block_heigth; r++)
				local_out[r] += A->bcsr_data[block*bs + c*A->bcsr_block_heigth + r] * vec_this_col;
		}
	}

	for(u64 r=0; r<A->bcsr_block_width; r++)
		y_data[target_block_row*A->bcsr_block_width + r] += local_out[r];
}

void SparseMatrixMultiplayer::sparseMVMultiplyCPU(SparseMatrixData* A, Vector* x, Vector* y)
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

		t[i] = std::thread(SparseMatrixMultiplayer::sparseMVMultiplyThreadHandler, &data[i]);
	}

	for(u64 i=0; i<n; i++)
		t[i].join();
}

i32 ceil_div(i32 x, i32 y)
{
	return 1 + ((x - 1) / y);
}

void SparseMatrixMultiplayer::sparseMMMultiplyGPU(SparseMatrixData* A, SparseMatrixData* B, SparseMatrixData* C)
{
	C->bcsr_col_idx.clear();
	C->bcsr_row_ptr.clear();

	C->bcsr_lines = A->bcsr_lines;
	C->bcsr_block_heigth = A->bcsr_block_heigth;
	C->bcsr_block_width  = A->bcsr_block_width;

	for(int i=0; i<C->bcsr_lines/C->bcsr_block_heigth + 1; i++)
		C->bcsr_row_ptr.push_back(0);

	HashMap* accumulator = new HashMap(A->bcsr_lines);

	std::vector<u64> teamRows = { 0, 1, 2, 3, 4, 5 };

	for(int phase=0; phase<2; phase++)
	{
		if(phase == 1)
		{
			C->bcsr_data.clear();
			C->bcsr_data.resize(C->bcsr_row_ptr[C->bcsr_row_ptr.size() - 1], 0.f);
			C->bcsr_col_idx.resize(C->bcsr_row_ptr[C->bcsr_row_ptr.size() - 1], 0);
		}

		for(u64 i : teamRows)
		{
			for(u64 a_idx = A->bcsr_row_ptr[i]; a_idx<A->bcsr_row_ptr[i+1]; a_idx++)
			{
				u64 j = A->bcsr_col_idx[a_idx];
				for(u64 b_idx = B->bcsr_row_ptr[j]; b_idx < B->bcsr_row_ptr[j+1]; b_idx++)
				{
					u64 col = B->bcsr_col_idx[b_idx];
					f32 tmp_val = B->bcsr_data[b_idx] * A->bcsr_data[a_idx];

					accumulator->insert(col, tmp_val);
				}
			}

			u64 len = accumulator->size();

			if(phase == 0)
			{
				C->bcsr_row_ptr[i+1] = C->bcsr_row_ptr[i] + len;
			}
		
			if(phase == 1)
			{
				sortKeyValuePair(accumulator->_keys, accumulator->_vals, len);

				for(u64 c_idx = C->bcsr_row_ptr[i]; c_idx < C->bcsr_row_ptr[i+1]; c_idx++)
				{
					C->bcsr_col_idx[c_idx]  = accumulator->_keys[c_idx - C->bcsr_row_ptr[i]];
					C->bcsr_data[c_idx] 		= accumulator->_vals[c_idx - C->bcsr_row_ptr[i]];
				}
			}

			accumulator->reset();
		}
	}

	delete accumulator;
}

void SparseMatrixMultiplayer::sparseMVMultiplyGPU(SparseMatrixData* A, Vector* x, Vector* y)
{
	compute::Buffer col_idx = compute::Buffer(A->bcsr_col_idx.data(), sizeof(u64)*A->bcsr_col_idx.size(), compute::Buffer::READ_WRITE, false);
	compute::Buffer row_ptr = compute::Buffer(A->bcsr_row_ptr.data(), sizeof(u64)*A->bcsr_row_ptr.size(), compute::Buffer::READ_WRITE, false);
	compute::Buffer A_buffe	= compute::Buffer(A->bcsr_data.data(), 		sizeof(f32)*A->bcsr_data.size(), 		compute::Buffer::READ_WRITE, false);
	compute::Buffer x_buffe = compute::Buffer(x->data(), 							sizeof(f32)*x->size(),	 						compute::Buffer::READ_WRITE, false);
	compute::Buffer y_buffe = compute::Buffer(sizeof(f32)*y->size(), compute::Buffer::READ_WRITE, compute::Buffer::MEM_GPU);

	// I think the same algorithm will work for the COO format
	// if the block dim is equal to 1x1=1, in that case the matrix
	// also need to have an row_ptr array with the same size
	// as the col_idx array and with the corresponding lines 
	const unsigned int block_dim = A->bcsr_lines/A->bcsr_block_heigth; //A->bcsr_block_heigth * A->bcsr_block_width;

	bsMV_kernel->setKernelArgument(0, sizeof(unsigned long int), &A->bcsr_block_heigth);
	bsMV_kernel->setKernelArgument(1, sizeof(unsigned long int), &A->bcsr_block_width);
	bsMV_kernel->setKernelArgument(2, BUFFER_ARG_SIZE, col_idx.upload());
	bsMV_kernel->setKernelArgument(3, BUFFER_ARG_SIZE, row_ptr.upload());
	bsMV_kernel->setKernelArgument(4, BUFFER_ARG_SIZE, A_buffe.upload());
	bsMV_kernel->setKernelArgument(5, BUFFER_ARG_SIZE, x_buffe.upload());
	bsMV_kernel->setKernelArgument(6, BUFFER_ARG_SIZE, y_buffe.upload());
	bsMV_kernel->setKernelArgument(7, 4*block_dim * sizeof(float), nullptr);
	bsMV_kernel->enqueue({ A->bcsr_data.size() }, { block_dim });

	f32* data = (f32*)y_buffe.download();

	for(u64 i=0; i<y->size(); i++)
	{
		y->data()[i] = data[i];
	}
}
