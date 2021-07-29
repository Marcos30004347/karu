#include <assert.h>

#include "algebra/linear/ConjugateGradient.hpp"

using namespace karu::algebra;

int conjugateGradient()
{
	Matrix A(2,2, {
		3, 2,
		2, 6,
	});

	Matrix x0(2, 1, {0,0});

	Matrix b(2, 1, {2,-8});

	Matrix x = conjGrad(A, x0, b, 0.000001);

	printf("solution:\n");
	printMatrix(x);

	Matrix B(2,2, {
		2, -1,
		-1, 2,
	});

	Matrix y0 = Matrix(2, 1, { 0, 0 });

	Matrix c(2, 1, { 1, 0 });

	Matrix y = conjGrad(B, y0, c, 0.0001);

	printf("solution:\n");
	printMatrix(y);

	Matrix C(5,5, {
		3, 1, 0, 0, 1,
		1, 4, 1, 1, 1,
		0, 1, 1, 1, 1,
		0, 1, 1, 1, 1,
		1, 1, 1, 1, 3,
	});

	Matrix z0 = Matrix(5, 1, { 0, 0, 0, 0, 0 });

	Matrix d(5, 1, { 5, 8, 4, 4, 7 });

	Matrix z = conjGrad(C, z0, d, 0.0000001);
	
	printf("solution:\n");	
	printMatrix(z);

	return 0;
}


int preconditionedConjugateGradient()
{
	Matrix I(2,2, {
		3, 0,
		0, 6,
	});

	Matrix A(2,2, {
		3, 2,
		2, 6,
	});

	Matrix x0(2, 1, {0, 0});

	Matrix b(2, 1, {2, -8});

	Matrix x = preConjGrad(A, x0, b, I, 0.000001);
	
	Matrix S = A*x;
	assert(S[0][0] == b[0][0]);
	assert(S[1][0] == b[1][0]);

	// printMatrix(S);

	// Matrix B(2,2, {
	// 	2, -1,
	// 	-1, 2,
	// });

	// Matrix y0 = Matrix(2, 1, { 0, 0 });

	// Matrix c(2, 1, { 1, 0 });

	// Matrix y = conjGrad(B, y0, c, 0.0001);

	// printf("solution:\n");
	// printMatrix(y);

	// Matrix C(5,5, {
	// 	3, 1, 0, 0, 1,
	// 	1, 4, 1, 1, 1,
	// 	0, 1, 1, 1, 1,
	// 	0, 1, 1, 1, 1,
	// 	1, 1, 1, 1, 3,
	// });

	// Matrix z0 = Matrix(5, 1, { 0, 0, 0, 0, 0 });

	// Matrix d(5, 1, { 5, 8, 4, 4, 7 });

	// Matrix z = conjGrad(C, z0, d, 0.0000001);
	
	// printf("solution:\n");	
	// printMatrix(z);

	return 0;
}

int main()
{
	preconditionedConjugateGradient();
}
