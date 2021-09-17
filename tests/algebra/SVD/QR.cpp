#include <assert.h>
#include "algebra/SVD/QR.hpp"

using namespace karu;
using namespace karu::algebra;

int main()
{
	Matrix A(4,4, {
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 1, 2, 3,
		4, 5, 6, 7,
	});
	Matrix Q;

	QR(A, Q, 2.22e-16);

	printMatrix(Q);
	printMatrix(A);
	printMatrix(Q*A);
	return 0;
}
