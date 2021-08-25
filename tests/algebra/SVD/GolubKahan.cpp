#include <assert.h>

#include "algebra/SVD/GolubKahan.hpp"
#include "algebra/SVD/SVD.hpp"

using namespace karu;
using namespace karu::algebra;

int main()
{
	// f32 diag[4] = {
	// 	1, 2, 3, 4,
	// };

	// f32 sup_diag[4] = {
	// 	0, 1, 1, 1,
	// };

	// Matrix uT = identity(5,5);
	// Matrix v = identity(4,4);

	// // golubKahanStep(diag, sup_diag, 5, 5, uT, v, 2.22045e-16);

	// golubKahanSVD(diag, sup_diag, 5, 4, uT, v, 2.22045e-16);

	Matrix A, U, V;
	f32 s[4];
	
	A = Matrix(4,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4
	});

	gsvd(A, U, s, V);
	
	A = Matrix(5,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4,
		1, 4, 7, 4
	});

	gsvd(A, U, s, V);

	A = Matrix(5,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4,
		1, 0, 0, 4
	});

	gsvd(A, U, s, V);


	return 0;
}
