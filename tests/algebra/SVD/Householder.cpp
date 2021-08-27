#include <assert.h>

#include "algebra/SVD/Householder.hpp"

using namespace karu;
using namespace karu::algebra;


void householderMatrixRowTest()
{
	i32 i;

	f32 beta;
	f32 x[4] = {1, 4, 7, 10};
	f32 v[4] = {0};
	f32 o[4] = {0};

	Matrix A(4, 3, {
		x[0], 2, 3,
		x[1], 5, 6,
		x[2], 8, 9,
		x[3], 11, 12
	});

	beta = houseCol(A, 4, 3, 0, 0, v, 2.22e-15);

	applyHouseholderToVector(v, beta, x, o, 4);
	
	assert(fabs(o[0] - norm(x, 4)) <= 2.22e-15);

	preHouseholderMatrix(v, beta, A, 4, 3, 0, 0);

	assert((fabs(A[0][0] - o[0])) <= 2.22e-15);
}

void householderMatrixColTest()
{
	i32 i;

	f32 beta;
	f32 x[3] = {1, 2, 3};
	f32 v[3] = {0};
	f32 o[3] = {0};

	Matrix A(4, 3, {
		x[0], x[1], x[2],
		4, 5, 6,
		5, 8, 9,
		6, 11, 12
	});

	beta = houseRow(A, 4, 3, 0, 0, v, 2.22e-15);

	applyHouseholderToVector(v, beta, x, o, 3);

	assert(fabs(o[0] - norm(x, 3)) <= 2.22e-15);

	posHouseholderMatrix(v, beta, A, 4, 3, 0, 0);

	assert((fabs(A[0][0] - o[0])) <= 2.22e-15);
}

void householderVectorTest()
{
	i32 i;

	const i32 n = 4;
	const i32 m = 3;

	f32 x[n] = {3, 2, 3, 4};
	f32 v[n] = {0};

	f32 Px[n] = {0};
	
	f32 beta = house(x, v, n, 2.22e-15);

	f32 out[n];

	applyHouseholderToVector(v, beta, x, out, n);

	assert((fabs(norm(x,n) - out[0])) < 2.22e-15);
}

void householderBidiagonalizationTests()
{
	Matrix A(4, 3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12
	});

	Matrix U, V;

	f32 diag[3];
	f32 sdiag[3];
	f32 trash;

	householderBidiagonalization(A, diag, sdiag, U, V, 2.22e-15);

	Matrix B(4, 3, {
		diag[0], sdiag[1], 0,
		0, diag[1], sdiag[2],
		0,       0,  diag[2],
		0,       0,        0
	});

	printMatrix(U);
	std::cout << "\n";
	printMatrix(B);
	std::cout << "\n";
	printMatrix(V);
	std::cout << "\n";
	printMatrix(U*B*transpose(V));
}

int main()
{
	householderVectorTest();
	householderMatrixRowTest();
	householderMatrixColTest();
	householderBidiagonalizationTests();
	return 0;
}
