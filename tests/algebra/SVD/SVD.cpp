#include <assert.h>

#include "algebra/SVD/Householder.hpp"
#include "algebra/SVD/GolubReinschSVD.hpp"
#include "algebra/SVD/SVD.hpp"
#include "algebra/SVD/GolubKahanLanczosBid.hpp"

using namespace karu;
using namespace karu::algebra;

void householderTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});

	Matrix u, v;

	f32 e[3];
	f32 q[3];

	golubReinschHouseholderBidiagonalization(a, q, e, u, v, 2.22e-16);

	std::cout << "e:\n";
	for(int i=0; i<3; i++)
		std::cout << e[i] << ", ";
	
	std::cout << "\n";
	std::cout << "q:\n";
	
	for(int i=0; i<3; i++)
		std::cout << q[i] << ", ";
	
	std::cout << "\n";
	std::cout << "u:\n";
	printMatrix(u);
	std::cout << "v:\n";
	printMatrix(v);
	std::cout << "a:\n";

	Matrix b(4, 3, {
		q[0], e[1], 0,
		0, q[1], e[2],
		0,   0,  q[2],
		0, 0, 0,
	});

	std::cout << "u*b*v':\n";
	printMatrix(u*b*transpose(v));
	std::cout << "diff':\n";
	printMatrix((u*b*transpose(v)) - a);

}


void bidiagonalTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});
	Matrix u, v;

	f32 gamma[3];
	f32 phi[3];
	f32 x;

	barlowBidiagonalization(a, gamma, phi, u, v, 2.22e-16, &x);

	for(i32 i=0; i<3; i++)
		std::cout << "gamma[i]: " << gamma[i] << "\n";

	for(i32 i=0; i<3; i++)
		std::cout << "phi[i]: " << phi[i] << "\n";

	printMatrix(u);

	Matrix B(3, 3, {
		gamma[0], phi[1], 0,
		0, gamma[1], phi[2],
		0, 		0, 	 gamma[2],
	});

	std::cout << "B:\n";
	printMatrix(B);
	std::cout << "U:\n";
	printMatrix(u);
	std::cout << "V:\n";
	printMatrix(v);

	std::cout << "U*B*V':\n";
	printMatrix(u*B*transpose(v));

	std::cout << "diff:\n";
	printMatrix(u*B*transpose(v) - Matrix(4,3, {1,2,3,4,5,6,7,8,9,10,11,12}));
}



void bidiagonalGolubKahanLaczosTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});
	
	f32 s[3], x;
	Matrix u, v;

	f32 eps = 1.e-15;
	f32 tol = 1.e-64/eps;

	f32 gamma[3];
	f32 phi[3];

	golubKahanLanczosBidiagonalization(a);
}

void bidiagonalGloubReinschSVDTests()
{
	Matrix a(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});
	
	f32 s[3], x;
	Matrix u, v;

	f32 eps = 1.e-15;
	f32 tol = 1.e-64/eps;

	f32 gamma[3];
	f32 phi[3];

	// To bidiagonal form
	golubReinschHouseholderBidiagonalization(a, gamma, phi, u, v, tol, &x);
	golubReinschBidiabonalSVD(gamma, phi, a.rows(), a.columns(), u, s, v, tol, eps, x);

	std::cout << "B\n";
	Matrix b(4, 3, {
		s[0], 0, 0,
		0, s[1], 0,
		0, 0,  s[2],
		0, 0,    0,
	});
	printMatrix(b);
	std::cout << "U\n";
	printMatrix(u);
	std::cout << "V'\n";
	printMatrix(transpose(v));
	std::cout << "\n";

	std::cout << "U*B*V':\n";
	printMatrix(u*b*transpose(v));
	std::cout << "residual:\n";
	printMatrix(u*b*transpose(v) - a);
	std::cout << "\n";

	barlowBidiagonalization(a, gamma, phi, u, v, tol, &x);
	printMatrix(a);
	
	golubReinschBidiabonalSVD(gamma, phi, a.rows(), a.columns(), u, s, v, tol, eps, x);

	b = Matrix(3, 3, {
		s[0], 0, 0,
		0, s[1], 0,
		0, 0,  s[2],
	});

	std::cout << "B\n";
	printMatrix(b);
	std::cout << "U\n";
	printMatrix(u);
	std::cout << "V'\n";
	printMatrix(v);
	std::cout << "\n";

	Matrix ubv = u*b*transpose(v);
	std::cout << "U*B*V':\n";
	printMatrix(ubv);
	std::cout << "residual:\n";
	printMatrix(ubv - a);

	f32 n = v[0][2]*v[0][2] + v[1][2]*v[1][2] + v[2][2]*v[2][2];
	n = sqrt(n);
	std::cout << n << "\n";

	f32 t = 1*v[0][2]*v[0][2] + 2*v[1][2]*v[1][2] + 3*v[2][2]*v[2][2];
	std::cout << t << "\n";

}

void svdTests()
{
	Matrix at, a, u, v, res;
	f32 s[3];

	a = Matrix(4,3, {
		1, 2, 3,
		4, 5, 6,
		7, 8, 9,
		10, 11, 12,
	});

	svd(a, u, s, v);
	res = u*diag(s, 4, 3)*transpose(v);

	for(i32 i=0; i<a.rows(); i++)
		for(i32 j=0; j<a.columns(); j++)
			assert(fabs(a[i][j] - res[i][j]) <= 2.22e-14);

	at = transpose(a);

	svd(at, u, s, v);

	res = u*transpose(diag(s, 4, 3))*transpose(v);

	for(i32 i=0; i<at.rows(); i++)
		for(i32 j=0; j<at.columns(); j++)
			assert(fabs(at[i][j] - res[i][j]) <= 2.22e-14);

}

int main()
{
	
	// householderTests();
	// bidiagonalGloubReinschSVDTests();
	// bidiagonalTests();
	bidiagonalGolubKahanLaczosTests();
	// svdTests();
	return 0;
}
