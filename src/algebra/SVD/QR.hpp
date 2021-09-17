#pragma once

#include "Householder.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"

using namespace karu;
using namespace karu::algebra;

// Algorithm 5.2.1 Householder QR from Matrix Computations 4th
void QR(Matrix& A, Matrix& Q, f32 tol)
{
	i32 i, j, k, q;
	
	i32 m = A.rows();
	i32 n = A.columns();

	f32* x = new f32[m];
	f32* v = new f32[m];
	f32* t = new f32[n];
	f32* b = new f32[n];

	for(j = 0; j < n; j++)
	{
		// x = A[j:m][j]
		for(i = j; i < m; i++)
		{
			x[i - j] = A[i][j];
		}

		b[j] = house(x, v, m - j, tol);

		// A[j:m][j:n] = (I - b*v*v') * A[j:m][j:n]
		// t[1:n-j] = v[1:j-m]' * A[j:m][j:n]
		for(i = 0; i < n - j; i++)
		{
			t[i] = 0;
			
			for(k = 0; k < m - j; k++)
			{
				t[i] += v[k] * A[j + k][j + i];
			}
		}
		// A[j:m][j:n] = A[j:m][j:n] -  b*v[1:j-m]*t[1:n-j]
		for(i = 0; i < m - j; i++)
		{
			for(k = 0; k < n - j; k++)
			{
				A[j + i][j + k] -= b[j] * v[i] * t[k];
			}
		}

		if(j < m - 1)
		{
			for(i = 0; i < m - (j + 1); i++)
			{
				A[j + 1 + i][j] = v[i + 1];
			}
		}
	}

	// Accumulate R matrix from Q matrix
	Q = identity(m, m);

	for(j = n - 1; j >= 0; j--)
	{
		v[j] = 1;
		for(i = 1; i < m - j; i++)
		{
			v[j + i] = A[j + i][j];
		}

		for(q = 0; q < m - j; q++)
		{
			t[j + q] = 0;
			for(i = 0; i < m - j; i++)
			{
				t[j + q] = t[j + q] + (v[j + i] * Q[j + i][j + q]);
			}
		}

		for(q = 0; q < m - j; q++)
		{
			for(i = 0; i < m - j; i++)
			{
				Q[j + q][j + i] -= b[j] * v[j + q] * t[j + i];
			}
		}
	}

	for(i = 1; i < m; i++)
	{
		for(j = 0; j < i; j++)
		{
			A[i][j] = 0.0;
		}
	}

	delete[] x;
	delete[] v;
	delete[] t;
	delete[] b;
}
