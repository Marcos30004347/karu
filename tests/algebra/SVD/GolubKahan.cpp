#include <assert.h>

#include "algebra/SVD/GolubKahan.hpp"

using namespace karu;
using namespace karu::algebra;

int main()
{
	f32 diag[4] = {
		1, 2, 3, 4
	};

	f32 sup_diag[4] = {
		0, 1, 1, 1
	};

	Matrix uT = identity(4,4);
	Matrix v = identity(4,4);

	// golubKahanStep(diag, sup_diag, 4, 4, uT, v, 2.22045e-16);

	golubKahanSVD(diag, sup_diag, 4, 4, uT, v, 2.22045e-16);

	return 0;
}
