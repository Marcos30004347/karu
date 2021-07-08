#include <assert.h>
#include <iostream>

#include "algebra/matrix/BlockSparseMatrixData.hpp"

using namespace karu;

int main()
{
	std::vector<f32> data;

	for(u64 i=0; i<36; i++)
	{
		data.push_back(i+1);
	}
	
	// M =
	// 1   2   3   0   0   0  7    8   9
	// 4   5   6   0   0   0  10  11  12
	// 0  13  14  15   0  19  20  21   0
	// 0  16  17  18   0  22  23  24   0
	// 0   0   0  25  26  27  31  32  33
	// 0   0   0  28  29  30  34  35  36
	
	algebra::BlockSparseMatrixData M = algebra::BlockSparseMatrixData(
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

	assert(M.get(0,0) == 1.f);
	assert(M.get(0,1) == 2.f);
	assert(M.get(0,2) == 3.f);
	assert(M.get(0,3) == 0.f);
	assert(M.get(0,4) == 0.f);
	assert(M.get(0,5) == 0.f);
	assert(M.get(0,6) == 7.f);
	assert(M.get(0,7) == 8.f);
	assert(M.get(0,8) == 9.f);

	assert(M.get(1,0) == 4.f);
	assert(M.get(1,1) == 5.f);
	assert(M.get(1,2) == 6.f);
	assert(M.get(1,3) == 0.f);
	assert(M.get(1,4) == 0.f);
	assert(M.get(1,5) == 0.f);
	assert(M.get(1,6) == 10.f);
	assert(M.get(1,7) == 11.f);
	assert(M.get(1,8) == 12.f);

	assert(M.get(2,0) == 0.f);
	assert(M.get(2,1) == 13.f);
	assert(M.get(2,2) == 14.f);
	assert(M.get(2,3) == 15.f);
	assert(M.get(2,4) == 0.f);
	assert(M.get(2,5) == 19.f);
	assert(M.get(2,6) == 20.f);
	assert(M.get(2,7) == 21.f);
	assert(M.get(2,8) == 0.f);

	assert(M.get(3,0) == 0.f);
	assert(M.get(3,1) == 16.f);
	assert(M.get(3,2) == 17.f);
	assert(M.get(3,3) == 18.f);
	assert(M.get(3,4) == 0.f);
	assert(M.get(3,5) == 22.f);
	assert(M.get(3,6) == 23.f);
	assert(M.get(3,7) == 24.f);
	assert(M.get(3,8) == 0.f);

	assert(M.get(4,0) == 0.f);
	assert(M.get(4,1) == 0.f);
	assert(M.get(4,2) == 0.f);
	assert(M.get(4,3) == 25.f);
	assert(M.get(4,4) == 26.f);
	assert(M.get(4,5) == 27.f);
	assert(M.get(4,6) == 31.f);
	assert(M.get(4,7) == 32.f);
	assert(M.get(4,8) == 33.f);

	assert(M.get(5,0) == 0.f);
	assert(M.get(5,1) == 0.f);
	assert(M.get(5,2) == 0.f);
	assert(M.get(5,3) == 28.f);
	assert(M.get(5,4) == 29.f);
	assert(M.get(5,5) == 30.f);
	assert(M.get(5,6) == 34.f);
	assert(M.get(5,7) == 35.f);
	assert(M.get(5,8) == 36.f);

	return 0;
}
