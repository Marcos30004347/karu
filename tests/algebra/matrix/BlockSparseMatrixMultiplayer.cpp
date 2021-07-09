#include <assert.h>
#include <iostream>

#include "algebra/matrix/BlockSparseMatrixData.hpp"
#include "algebra/matrix/BlockSparseMatrixMultiplayer.hpp"
#include "algebra/core/compute/Compute.hpp"

using namespace karu;

int main()
{
	algebra::compute::Context::initContext();

	std::vector<f32> data;

	for(u64 i=0; i<36; i++)
	{
		data.push_back(i+1);
	}
	
	algebra::BlockSparseMatrixData M = algebra::BlockSparseMatrixData(
		3, 2, // blocks will have 3 elements width, 2 elements heigth
		6, 9, // 6 lines and 9 columns
		{0, 2, 4, 6}, // row_ptr
		{0, 6, 1, 5, 3, 6}, // col_idx
		{
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f,
			18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f,
			26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f,
			34.f, 35.f, 36.f
		}
	);

	algebra::BlockSparseMatrixData x = algebra::BlockSparseMatrixData(
		1, 2, // blocks will have 1 elements width, 6 elements heigth
		9, 1, // 6 lines and 9 columns
		{0, 1, 2, 3, 4, 5}, // row_ptr
		{0, 0, 0, 0}, // col_idx
		{
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 0.f
		}
	);

	algebra::BlockSparseMatrixData y = algebra::BlockSparseMatrixData(
		1, 2, // blocks will have 1 elements width, 6 elements heigth
		6, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{
			0.f, 0.f, 0.f, 0.f, 0.f,
		}
	);

	algebra::BlockSparseMatrixMultiplayer::sparseMVMultiplyThreaded(&M, &x, &y);

	assert(y.get(0, 0) == 242.f);
	assert(y.get(1, 0) == 272.f);
	assert(y.get(2, 0) == 584.f);
	assert(y.get(3, 0) == 614.f);
	assert(y.get(4, 0) == 1205.f);
	assert(y.get(5, 0) == 1244.f);

	algebra::BlockSparseMatrixData A = algebra::BlockSparseMatrixData(
		2, 2, // blocks will have 3 elements width, 2 elements heigth
		6, 6, // 6 lines and 9 columns
		{0, 2, 3, 5}, // row_ptr
		{0, 4, 2, 0, 4}, // col_idx
		{
			1.f, 2.f, 3.f, 4.f,
			5.f, 6.f, 7.f, 8.f,
			9.f, 10.f, 11.f, 12.f,
			13.f, 14.f, 15.f, 16.f,
			17.f, 18.f, 19.f, 20.f,
		}
	);

	algebra::BlockSparseMatrixData z = algebra::BlockSparseMatrixData(
		1, 2, // blocks will have 1 elements width, 6 elements heigth
		6, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f}
	);

	algebra::BlockSparseMatrixData w = algebra::BlockSparseMatrixData(
		1, 2, // blocks will have 1 elements width, 6 elements heigth
		6, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}
	);

	algebra::BlockSparseMatrixMultiplayer::sparseMVMultiplyThreaded(&A, &z, &w);

	assert(w.get(0, 0) == 74.f);
	assert(w.get(1, 0) == 88.f);
	assert(w.get(2, 0) == 71.f);
	assert(w.get(3, 0) == 78.f);
	assert(w.get(4, 0) == 242.f);
	assert(w.get(5, 0) == 256.f);

	algebra::BlockSparseMatrixData B = algebra::BlockSparseMatrixData(
		3, 3, // blocks will have 3 elements width, 2 elements heigth
		9, 9, // 6 lines and 9 columns
		{0, 3, 6, 9 }, // row_ptr
		{0, 3, 6, 0, 3, 6, 0, 3, 6}, // col_idx
		{
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
		}
	);

	algebra::BlockSparseMatrixData o = algebra::BlockSparseMatrixData(
		1, 3, // blocks will have 1 elements width, 6 elements heigth
		9, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f }
	);

	algebra::BlockSparseMatrixData u = algebra::BlockSparseMatrixData(
		1, 3, // blocks will have 1 elements width, 6 elements heigth
		9, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }
	);

	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<9; j++)
		{
			std::cout << B.get(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<1; j++)
		{
			std::cout << o.get(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	algebra::BlockSparseMatrixMultiplayer::sparseMVMultiplyGPU(&B, &o, &u);

	std::cout << std::endl;
	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<1; j++)
		{
			std::cout << u.get(i, j) << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	algebra::BlockSparseMatrixData C = algebra::BlockSparseMatrixData(
		3, 3, // blocks will have 3 elements width, 2 elements heigth
		9, 9, // 6 lines and 9 columns
		{0, 2, 3, 5 }, // row_ptr
		{0, 6, 3, 0, 6}, // col_idx
		{
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
			1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
		}
	);

	algebra::BlockSparseMatrixData q = algebra::BlockSparseMatrixData(
		1, 3, // blocks will have 1 elements width, 6 elements heigth
		9, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f }
	);

	algebra::BlockSparseMatrixData t = algebra::BlockSparseMatrixData(
		1, 3, // blocks will have 1 elements width, 6 elements heigth
		9, 1, // 6 lines and 9 columns
		{0, 1, 2, 3}, // row_ptr
		{0, 0, 0}, // col_idx
		{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }
	);

	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<9; j++)
		{
			std::cout << C.get(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;
	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<1; j++)
		{
			std::cout << q.get(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	algebra::BlockSparseMatrixMultiplayer::sparseMVMultiplyGPU(&C, &q, &t);

	for(u64 i=0; i<9; i++)
	{
		for(u64 j=0; j<1; j++)
		{
			std::cout << t.get(i, j) << " ";
		}
		std::cout << std::endl;
	}

	std::cout << std::endl;

	algebra::compute::Context::stopContext();
	return 0;
}

