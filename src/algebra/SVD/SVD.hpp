#pragma once


#include "algebra/SVD/GolubReinschSVD.hpp"
#include "algebra/SVD/Householder.hpp"
#include "algebra/SVD/GolubKahan.hpp"

namespace karu::algebra {

void reorder(f32* s, u32 n, Matrix* matrices[2])
{
	for (int i = 0; i < n; ++i) 
	{
		double s_last = s[i];
		int i_last = i;
	
		for (int j = i + 1; j < n; ++j)
		{
			if (s[j] >  s_last)
			{
				s_last = s[j];
				i_last = j;
			}
		}
		
		if (i_last != i)
		{
			double tmp;
			tmp = s[i];
			s[i] = s[i_last];
			s[i_last] = tmp;

			for (int j = 0; j < 2; ++j)
			{
				int rows = matrices[j]->rows();
				int cols = matrices[j]->columns();
			
				Matrix* M = matrices[j];
			
				for (int k = 0; k < rows; ++k) 
				{
					tmp = M->m_data.get(k, i);
					M->m_data.set(k, i, M->m_data.get(k, i_last));
					M->m_data.set(k, i_last, tmp);
				}
			}
		}
	}
}

void svd(Matrix& a, Matrix& u, f32* s, Matrix& v)
{
	f32 x, eps = 1.e-15, tol = 1.e-64/eps;

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

	f32* gamma = new f32[n];
	f32* phi = new f32[n];

	golubReinschHouseholderBidiagonalization(A, gamma, phi, u, v, tol, &x);

	for(i32 i=0; i<n; i++)
		std::cout << gamma[i] << " ";
	std::cout << "\n";
	for(i32 i=0; i<n; i++)
		std::cout << phi[i] << " ";
	std::cout << "\n";

	golubReinschBidiabonalSVD(gamma, phi, m, n, u, s, v, tol, eps, x);

	if(A.rows() == a.columns() && A.columns() == a.rows())
	{
		Matrix tmp = u;
		u = v;
		v = tmp;
	}

	delete phi;
	delete gamma;

	Matrix* vecs[2];
	vecs[0] = &u;
	vecs[1] = &v;
	reorder(s, n, vecs);
}


void gsvd(Matrix& a, Matrix& u, f32* s, Matrix& v, f32 tol = 2.22e-15)
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

	f32* gamma = new f32[n];

	f32* phi   = new f32[n];

	phi[0] = 0;

	householderBidiagonalization(A, gamma, phi, u, v, tol);

	printMatrix(u);
	std::cout << "\n";
	printMatrix(v);
	std::cout << "\n";
	
	for(i32 i=0; i<n; i++)
		std::cout << gamma[i] << " ";
	std::cout << "\n";

	for(i32 i=0; i<n; i++)
		std::cout << phi[i] << " ";
	std::cout << "\n";

	Matrix B(4,4, {
		gamma[0], phi[1], 0, 0,
		0, gamma[1], phi[2], 0,
		0, 0, gamma[2], phi[3],
		0, 0, 0, gamma[3]
	});

	printMatrix(u * B * transpose(v));
	// golubKahanSVD(gamma, phi, m, n, u, v, tol);

	// for(i32 i=0; i<n; i++)
	// 	s[i] = gamma[i];

	// if(A.rows() == a.columns() && A.columns() == a.rows())
	// {
	// 	Matrix& tmp = u;
	// 	u = v;
	// 	v = tmp;
	// }

	delete[] phi;
	delete[] gamma;
}


}
