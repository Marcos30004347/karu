#pragma once


#include "algebra/SVD/GolubReinschSVD.hpp"
#include "algebra/SVD/Householder.hpp"
#include "algebra/SVD/GolubKahan.hpp"

namespace karu::algebra {


void svd(Matrix& a, Matrix& u, f32* s, Matrix& v, f32 tol = 2.22e-15)
{
	f32 x, eps = 2.22e-16;

	Matrix A;


	if(a.columns() > a.rows())
	{
		A = transpose(a);
	}
	else
	{
		A = Matrix(a);
	}

	u32 m = A.rows();
	u32 n = A.columns();

	u = identity(m,m);
	v = identity(n,n);


	f32* phi   = new f32[n];
	phi[0] = 0;

	printMatrix(A, 5);

	householderBidiagonalization(A, s, phi, u, v, tol);
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";
	std::cout << "\n";

	printMatrix(bidiag(s, phi, m, n), 5, 2.2e-5);

	printMatrix(u*bidiag(s, phi, m, n)*transpose(v), 5, 2.2e-5);

	std::cout << "\n";

	golubKahanSVD(s, phi, m, n, u, v, tol);

	// printMatrix(u*diag(s, m, n)*transpose(v), 5, 2.2e-5);
	// std::cout << "\n";

	if(a.columns() != a.rows() && A.rows() == a.columns() && A.columns() == a.rows())
	{
		Matrix tmp = u;
		u = v;
		v = tmp;
	}


	delete[] phi;
}


}
