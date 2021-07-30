#include <assert.h>
#include <iostream>

#include "algebra/sparse/SparseMatrixData.hpp"
#include "algebra/sparse/SparseMatrixMultiplayer.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;

int main()
{
	algebra::compute::Context::initContext();

	// algebra::SparseMatrixData M = algebra::SparseMatrixData(
	// 	3, 2, // blocks will have 3 elements width, 2 elements heigth
	// 	6, 9, // 6 lines and 9 columns
	// 	{0, 2, 4, 6}, // row_ptr
	// 	{0, 6, 1, 5, 3, 6}, // col_idx
	// 	{
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		10.f, 11.f, 12.f, 13.f, 14.f, 15.f, 16.f, 17.f,
	// 		18.f, 19.f, 20.f, 21.f, 22.f, 23.f, 24.f, 25.f,
	// 		26.f, 27.f, 28.f, 29.f, 30.f, 31.f, 32.f, 33.f,
	// 		34.f, 35.f, 36.f
	// 	}
	// );

	// algebra::Vector x = algebra::Vector({1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 0.f});
	// algebra::Vector y = algebra::Vector({0.f, 0.f, 0.f, 0.f, 0.f});
	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyCPU(&M, &x, &y);

	// assert(y[0] == 242.f);
	// assert(y[1] == 272.f);
	// assert(y[2] == 584.f);
	// assert(y[3] == 614.f);
	// assert(y[4] == 1205.f);
	// assert(y[5] == 1244.f);

	// algebra::SparseMatrixData M0 = algebra::SparseMatrixData(3,3, 9,9);

	// M0.defineMatrixSize(9,9);
	// M0.defineBlockSize(3,3);

	// M0.defineElement(0, 1);

	// algebra::SparseMatrixData M0 = algebra::SparseMatrixData(
	// 	3, 3, // blocks will have 3 elements width, 2 elements heigth
	// 	9, 9, // 6 lines and 9 columns
	// 	{0, 3, 6, 9 }, // row_ptr
	// 	{0, 3, 6, 0, 3, 6, 0, 3, 6}, // col_idx
	// 	{
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 	}
	// );

	// // algebra::Vector v0 = algebra::Vector({ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f });
	// // algebra::Vector v1 = algebra::Vector({ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f });

	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyCPU(&M0, &v0, &v1);
	// assert(v1 == algebra::Vector({198, 243, 288,198,243,288,198,243,288}));
	
	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyGPU(&M0, &v0, &v1);
	// assert(v1 == algebra::Vector({198, 243, 288,198,243,288,198,243,288}));

	// algebra::SparseMatrixData M1 = algebra::SparseMatrixData(
	// 	1, 1,
	// 	9, 9,
	// 	{ 0, 9, 18, 27, 36, 45, 54, 63, 72, 81 }, // row_ptr
	// 	{
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8,
	// 		0,1,2,3,4,5,6,7,8
	// 	}, // col_idx
	// 	{
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 	}
	// );

	// algebra::Vector v2 = algebra::Vector({ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f });
	// algebra::Vector v3 = algebra::Vector({ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f });

	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyCPU(&M1, &v2, &v3);
	// assert(v3 == algebra::Vector({ 285, 285, 285, 285, 285, 285, 285, 285,285 }));

	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyGPU(&M1, &v2, &v3);
	// assert(v3 == algebra::Vector({ 285, 285, 285, 285, 285, 285, 285, 285,285 }));

	// algebra::SparseMatrixData A = algebra::SparseMatrixData(
	// 	2, 2, // blocks will have 3 elements width, 2 elements heigth
	// 	6, 6, // 6 lines and 9 columns
	// 	{0, 2, 3, 5}, // row_ptr
	// 	{0, 4, 2, 0, 4}, // col_idx
	// 	{
	// 		1.f, 2.f, 3.f, 4.f,
	// 		5.f, 6.f, 7.f, 8.f,
	// 		9.f, 10.f, 11.f, 12.f,
	// 		13.f, 14.f, 15.f, 16.f,
	// 		17.f, 18.f, 19.f, 20.f,
	// 	}
	// );

	// algebra::SparseMatrixData z = algebra::SparseMatrixData(
	// 	1, 2, // blocks will have 1 elements width, 6 elements heigth
	// 	6, 1, // 6 lines and 9 columns
	// 	{0, 1, 2, 3}, // row_ptr
	// 	{0, 0, 0}, // col_idx
	// 	{ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f}
	// );

	// algebra::SparseMatrixData w = algebra::SparseMatrixData(
	// 	1, 2, // blocks will have 1 elements width, 6 elements heigth
	// 	6, 1, // 6 lines and 9 columns
	// 	{0, 1, 2, 3}, // row_ptr
	// 	{0, 0, 0}, // col_idx
	// 	{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}
	// );

	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyThreaded(&A, &z, &w);

	// assert(w.get(0, 0) == 74.f);
	// assert(w.get(1, 0) == 88.f);
	// assert(w.get(2, 0) == 71.f);
	// assert(w.get(3, 0) == 78.f);
	// assert(w.get(4, 0) == 242.f);
	// assert(w.get(5, 0) == 256.f);



	// algebra::SparseMatrixData C = algebra::SparseMatrixData(
	// 	3, 3, // blocks will have 3 elements width, 2 elements heigth
	// 	9, 9, // 6 lines and 9 columns
	// 	{0, 2, 3, 5 }, // row_ptr
	// 	{0, 6, 3, 0, 6}, // col_idx
	// 	{
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 		1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f,
	// 	}
	// );

	// algebra::SparseMatrixData q = algebra::SparseMatrixData(
	// 	1, 3, // blocks will have 1 elements width, 6 elements heigth
	// 	9, 1, // 6 lines and 9 columns
	// 	{0, 1, 2, 3}, // row_ptr
	// 	{0, 0, 0}, // col_idx
	// 	{ 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f }
	// );

	// algebra::SparseMatrixData t = algebra::SparseMatrixData(
	// 	1, 3, // blocks will have 1 elements width, 6 elements heigth
	// 	9, 1, // 6 lines and 9 columns
	// 	{0, 1, 2, 3}, // row_ptr
	// 	{0, 0, 0}, // col_idx
	// 	{ 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }
	// );


	// algebra::SparseMatrixMultiplayer::sparseMVMultiplyGPU(&C, &q, &t);

	// algebra::SparseMatrixData A0 = algebra::SparseMatrixData(
	// 	1, 1,
	// 	6, 6,
	// 	{ 0, 2, 4, 6, 9, 10, 11 }, // row_ptr
	// 	{ 0, 5, 3, 4, 1, 2, 1, 2, 3, 0, 5 }, // col_idx
	// 	{ 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f }
	// );

	// A0.print();

	// std::cout << std::endl;

	// algebra::SparseMatrixData B0 = algebra::SparseMatrixData(
	// 	1, 1, // blocks will have 3 elements width, 2 elements heigth
	// 	6, 6, // 6 lines and 9 columns
	// 	{ 0, 2, 4, 6, 8, 10, 16 }, // row_ptr
	// 	{ 0, 5, 2, 5, 1, 3, 3, 4, 4, 5, 0, 1, 2, 3, 4, 5 }, // col_idx
	// 	{ 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f }
	// );

	// B0.print();
	// std::cout << std::endl;
	
	// algebra::SparseMatrixData C0 = algebra::SparseMatrixData(
	// 	1, 1,
	// 	6, 6,
	// 	{},
	// 	{},
	// 	{}
	// );

	// C0.print();
	// std::cout << std::endl;

	algebra::compute::Context::stopContext();

	return 0;
}

