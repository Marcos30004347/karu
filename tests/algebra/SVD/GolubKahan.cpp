#include <assert.h>

#include "algebra/SVD/GolubKahan.hpp"
#include "algebra/SVD/SVD.hpp"

using namespace karu;
using namespace karu::algebra;

void assertMatrixIsClose(Matrix A, Matrix B, f32 tol = 2.22e-14)
{
	for(i32 i=0; i<A.rows(); i++)
	{
		for(i32 j=0; j<A.columns(); j++)
		{
			assert(fabs(A[i][j] - B[i][j]) <= tol);
		}
	}
}

int main()
{
	Matrix A, U, V;

	f32* singular_values = new f32[4];

	A = Matrix(4,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values, A.rows(), A.columns())*V);

	A = Matrix(5,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4,
		1, 4, 7, 4
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values, A.rows(), A.columns())*V);

	A = Matrix(5,4, {
		1, 1, 0, 0,
		0, 2, 1, 0,
		0, 2, 3, 1,
		1, 0, 0, 4,
		1, 0, 0, 4
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values, A.rows(), A.columns())*V);

	A = Matrix(5,4, {
		1, 1, 3, 7,
		4, 2, 9, 4,
		3, 2, 3, 1,
		1, 4, 1, 4,
		1, 5, 3, 4
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values, A.rows(), A.columns())*V);

	delete[] singular_values;

	singular_values = new f32[5];
	
	A = Matrix(4,5, {
		1, 1, 3, 7, 2,
		4, 2, 9, 4, 3,
		3, 2, 3, 1, 7,
		1, 4, 1, 4, 9,
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values, A.rows(), A.columns())*V);

	delete[] singular_values;

	singular_values = new f32[9];

	A = Matrix(9,9, {
    1.10, 	1.15, 	33.4, 	7.19, 	2.31, 	3.567, 	4.67, 	10.4, 	1,
    4.67, 	6.75, 	22.7, 	1.89, 	4.21, 	9.537, 	9.67, 	10.4, 	1,
    9.47, 	4.55, 	13.4, 	8.56, 	7.21, 	10.537, 14.67, 	703.4, 	1,
    17.37, 	8.51, 	53.4, 	77.56, 	9.35, 	15.517, 19.37, 	713.4, 	1,
    17.37, 	8.51, 	53.4, 	77.56, 	9.35, 	15.517, 19.37, 	713.4, 	1,
    2.55, 	18.51, 	52.5, 	574.53, 0.35, 	11.17, 	1.47, 	73.4, 	1,
    33.55, 	21.51, 	15.5, 	513.55, 9.35, 	14.17, 	77.45, 	94.1, 	1,
    18.15, 	34.11, 	24.5, 	51.55, 	16.35, 	13.17, 	123.45, 18.991, 1,
    18.15, 	34.11, 	24.5, 	51.55, 	16.35, 	13.17, 	123.45, 18.991, 1,
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values,  A.rows(), A.columns())*V,  2.22e-12);

	A = Matrix(8,9, {
    1.10, 	1.15, 	33.4, 	7.19, 	2.31, 	3.567, 	4.67, 	10.4, 	1,
    4.67, 	6.75, 	22.7, 	1.89, 	4.21, 	9.537, 	9.67, 	10.4, 	1,
    9.47, 	4.55, 	13.4, 	8.56, 	7.21, 	10.537, 14.67, 	703.4, 	1,
    17.37, 	8.51, 	53.4, 	77.56, 	9.35, 	15.517, 19.37, 	713.4, 	1,
    17.37, 	8.51, 	53.4, 	77.56, 	9.35, 	15.517, 19.37, 	713.4, 	1,
    2.55, 	18.51, 	52.5, 	574.53, 0.35, 	11.17, 	1.47, 	73.4, 	1,
    33.55, 	21.51, 	15.5, 	513.55, 9.35, 	14.17, 	77.45, 	94.1, 	1,
    18.15, 	34.11, 	24.5, 	51.55, 	16.35, 	13.17, 	123.45, 18.991, 1,
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values,  A.rows(), A.columns())*V,  2.22e-12);

	delete[] singular_values;

	singular_values = new f32[3];

	A = Matrix(3,3, {
		1.2208058232190191e-07,   5.1200019258368155e+00,   2.5927050417706710e-07, 
			5.1200018037558994e+00,  -3.8135154647843931e-07,   5.1200018037561046e+00, 
		-1.2208057932616700e-07,  -5.1200016816754985e+00,  -2.5927049462914908e-07,
	});

	svd(A, U, singular_values, V);

	assertMatrixIsClose(A, U*diag(singular_values,  A.rows(), A.columns())*V,  2.22e-12);
	
	delete[] singular_values;

	return 0;
}
