#include <assert.h>

#include "algebra/SVD/GolubKahan.hpp"

using namespace karu;
using namespace karu::algebra;

int main()
{
	f32 diag[5] = {
		1, 2, 3, 4, 7,
	};

	f32 sup_diag[5] = {
		0, 1, 1, 1, 2,
	};

	Matrix uT = identity(5,5);
	Matrix v = identity(5,5);

	// golubKahanStep(diag, sup_diag, 5, 5, uT, v, 2.22045e-16);

	golubKahanSVD(diag, sup_diag, 5, 5, uT, v, 2.22045e-16);

	return 0;
}
