#include <assert.h>
#include <iostream>

#include "algebra/sparse/SparseMatrixData.hpp"

using namespace karu;

int main()
{
	std::vector<f32> data;

	for(u64 i=0; i<36; i++)
	{
		data.push_back(i+1);
	}
	
	algebra::SparseMatrixData M = algebra::SparseMatrixData(
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

	// for(u64 i=0; i<6; i++)
	// {
	// 	for(u64 j=0; j<9; j++)
	// 	{
	// 		std::cout << M.get(i, j) << " ";
	// 	}
	// 	std::cout << std::endl;
	// }

	assert(M.get(0,0) == 1.f);
	assert(M.get(0,1) == 3.f);
	assert(M.get(0,2) == 5.f);
	assert(M.get(0,3) == 0.f);
	assert(M.get(0,4) == 0.f);
	assert(M.get(0,5) == 0.f);
	assert(M.get(0,6) == 7.f);
	assert(M.get(0,7) == 9.f);
	assert(M.get(0,8) == 11.f);

	assert(M.get(1,0) == 2.f);
	assert(M.get(1,1) == 4.f);
	assert(M.get(1,2) == 6.f);
	assert(M.get(1,3) == 0.f);
	assert(M.get(1,4) == 0.f);
	assert(M.get(1,5) == 0.f);
	assert(M.get(1,6) == 8.f);
	assert(M.get(1,7) == 10.f);
	assert(M.get(1,8) == 12.f);

	assert(M.get(2,0) == 0.f);
	assert(M.get(2,1) == 13.f);
	assert(M.get(2,2) == 15.f);
	assert(M.get(2,3) == 17.f);
	assert(M.get(2,4) == 0.f);
	assert(M.get(2,5) == 19.f);
	assert(M.get(2,6) == 21.f);
	assert(M.get(2,7) == 23.f);
	assert(M.get(2,8) == 0.f);

	assert(M.get(3,0) == 0.f);
	assert(M.get(3,1) == 14.f);
	assert(M.get(3,2) == 16.f);
	assert(M.get(3,3) == 18.f);
	assert(M.get(3,4) == 0.f);
	assert(M.get(3,5) == 20.f);
	assert(M.get(3,6) == 22.f);
	assert(M.get(3,7) == 24.f);
	assert(M.get(3,8) == 0.f);

	assert(M.get(4,0) == 0.f);
	assert(M.get(4,1) == 0.f);
	assert(M.get(4,2) == 0.f);
	assert(M.get(4,3) == 25.f);
	assert(M.get(4,4) == 27.f);
	assert(M.get(4,5) == 29.f);
	assert(M.get(4,6) == 31.f);
	assert(M.get(4,7) == 33.f);
	assert(M.get(4,8) == 35.f);

	assert(M.get(5,0) == 0.f);
	assert(M.get(5,1) == 0.f);
	assert(M.get(5,2) == 0.f);
	assert(M.get(5,3) == 26.f);
	assert(M.get(5,4) == 28.f);
	assert(M.get(5,5) == 30.f);
	assert(M.get(5,6) == 32.f);
	assert(M.get(5,7) == 34.f);
	assert(M.get(5,8) == 36.f);

	return 0;
}
