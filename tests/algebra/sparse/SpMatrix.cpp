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
 
	printMatrix(A);
	printMatrix(C);
	
	algebra::compute::Context::stopContext();
	return 0;
}
