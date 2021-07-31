#include <assert.h>
#include "algebra/sparse/SpMatrix.hpp"
#include "algebra/compute/Compute.hpp"

using namespace karu;
using namespace karu::algebra;

int main()
{
	algebra::compute::Context::initContext();

	// 1 2 3
	// 4 5 6
	// 7 8 9
	SpMatrix A(
		9,9, 3,3,
		{
			0, 2, 3, 5 // rows ptr
		},
		{
			0, 6, 3, 0, 6 // columns idx
		},
		{
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block0
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block1
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block2
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block3
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block4
		}
	);

	Matrix B(9, 1, { 1, 2, 3, 4, 5, 6, 7, 8, 9 });

	Matrix C = A*B;

	assert(C[0][0] == 64);
	assert(C[1][0] == 154);
	assert(C[2][0] == 244);
	assert(C[3][0] == 32);
	assert(C[4][0] == 77);
	assert(C[5][0] == 122);
	assert(C[6][0] == 64);
	assert(C[7][0] == 154);
	assert(C[8][0] == 244);

	SpMatrix D(
		9,9, 3,3,
		{
			0, 2, 3, 5 // rows ptr
		},
		{
			0, 6, 3, 0, 6 // columns idx
		},
		{
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block0
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block1
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block2
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block3
			1, 4, 7, 2, 5, 8, 3, 6, 9, // block4
		}
	);

	Matrix E(9, 1, { 1, 2, 3, 4, 5, 6, 7, 8, 9 }, 3, 3);

	Matrix F = A*B;

	assert(F[0][0] == 64);
	assert(F[1][0] == 154);
	assert(F[2][0] == 244);
	assert(F[3][0] == 32);
	assert(F[4][0] == 77);
	assert(F[5][0] == 122);
	assert(F[6][0] == 64);
	assert(F[7][0] == 154);
	assert(F[8][0] == 244);

	algebra::compute::Context::stopContext();
	
	return 0;
}
