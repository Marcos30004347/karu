#include "Householder.hpp"
#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"

using namespace karu;
using namespace karu::algebra;

// Algorithm 7.4.2 from Matrix Computations 4th
void hessenbergForm(Matrix& A, f32 tol)
{
	i32 k, i, j, n;
	f32 b, *v, *t, *x;

	n = A.rows();

	v = new f32[n];
	t = new f32[n];
	x = new f32[n];

	for(k = 0; k < n - 2; k++)
	{
		for(j = k + 1; j < n; j++)
		{
			x[j - (k + 1)] = A[j][k];
		}

		b = house(x, v, n - (k + 1), tol);

		// A[k + 1 : n, k:n] = (I - b * v*v')*A[k + 1 : n, k:n]
		// t[1:n - k] = v[1:n - (k + 1)]' * A[k + 1 : n, k:n]
		for(i = 0; i < n - k; i++)
		{
			t[i] = 0;
			for(j = 0; j < n - (k + 1); j++)
			{
				t[i] += v[j] * A[k + 1 + j][i + k];
			}
		}
		//	A[k + 1 : n, k:n]	= A[k + 1 : n, k:n] - b * v[1:n - (k + 1)] * t[1:n - k]
		for(i = 0; i < n - k; i++)
		{
			for(j = 0; j < n - (k + 1); j++)
			{
				A[k + 1 + j][i + k] = A[k + 1 + j][i + k] - b * v[i] * t[j];
			}
		}
		// Store esential part of u where new zeros where introduced.
		for(i = 0; i < n - (k + 1); i++)
		{
			A[i + (k + 2)][k] = v[i + 1];
		}

		// A[1:n, k + 1:n] = A[1:n, k + 1:n] * (I - b * v*v')
		// t[1:n] = A[1:n, k + 1:n] * v[1:n - (k + 1)]
		for(i = 0; i < n; i++)
		{
			t[i] = 0;
			for(j = 0; j < n - (k + 1); j++)
			{
				t[i] += v[j] * A[i][k + 1 + j];
			}
		}
		//	A[1:n, k + 1:n]	= A[1:n, k + 1:n] - b * t[1:n] * v[1:n - (k + 1)]
		for(i = 0; i < n; i++)
		{
			for(j = 0; j < n - (k + 1); j++)
			{
				A[i][k + 1 + j] = A[i][k + 1 + j] - b * t[i] * v[j];
			}
		}
		// Store esential part of u where new zeros where introduced.
		for(i = 0; i < n - (k + 2); i++)
		{
			A[0][i + (k + 2)] = v[i + 1];
		}
	}

	delete[] v;
	delete[] t;
	delete[] x;
}

void francisQRStep(Matrix& H, f32 tol)
{
	i32 m, n, k, i, j;
	f32 b, q, r, t, s, x, y, z, u[3], v[3], *l;

	n = H.rows();
	m = n - 1;

	l = new f32[n];

	s = H[m - 1][m - 1] + H[n - 1][n - 1];
	t = H[m - 1][m - 1] - H[m - 1][n - 1] * H[n - 1][m - 1];

	x = H[0][0] * H[0][0] + H[0][1] * H[1][0] - s * H[0][0] + t;
	y = H[1][0] * (H[0][0] + H[1][1] - s);
	z = H[1][0] * H[2][1];

	for(k = 0; k < n - 3; k++)
	{
		u[0] = x; 
		u[1] = y; 
		u[2] = z;
	
		b = house(u, v, 3, tol);

		q = std::max(1, k);

		// H[k + 1 : k + 3, q:n] = (I - b * v*v')*H[k + 1 : k + 3, q:n]
		// l[1:n - q] = v[1:3]' * H[k + 1 : k + 3, q:n]
		for(i = 0; i < n - q; i++)
		{
			l[i]  = v[0] * H[k + 0][i + q];
			l[i] += v[1] * H[k + 1][i + q]; 
			l[i] += v[2] * H[k + 2][i + q];
		}
		// H[k + 1 : k + 3, q:n] = H[k + 1:k + 3][q:n] - b * v[1:3] * l[1:n - q]
		for(i = 0; i < n - q; i++)
		{
			H[k + 0][i + q] = H[k + 0][i + q] - b * v[0] * l[i];
			H[k + 1][i + q] = H[k + 1][i + q] - b * v[1] * l[i];
			H[k + 2][i + q] = H[k + 2][i + q] - b * v[2] * l[i];
		}
		// Store esential part of u where new zeros where introduced.
		for(i = 0; i < n - (q + 1); i++)
		{
			H[i + k + 1][q] = v[i + 1];
		}

		r = std::min(k + 4, n);
	
		// H[1 : r, k + 1 : k + 3] = H[1 : r, k + 1 : k + 3] * (I - b * v*v')
		// l[1:r] = H[1 : r, k + 1 : k + 3] * v[1:3]
		for(i = 0; i < r; i++)
		{
			l[i]  = H[i][k + 0] * v[0];
			l[i] += H[i][k + 1] * v[1]; 
			l[i] += H[i][k + 2] * v[2]; 
		}
		// H[1 : r, k + 1 : k + 3] = H[1 : r, k + 1 : k + 3] - b * l[1:r] * v[1:3]'
		for(i = 0; i < r; i++)
		{
			H[i][k + 0] = H[i][k + 0] - b * l[i] * v[0];
			H[i][k + 1] = H[i][k + 1] - b * l[i] * v[1];
			H[i][k + 2] = H[i][k + 2] - b * l[i] * v[2];
		}
		// Store esential part of u where new zeros where introduced.
		H[0][k + 1] = v[1];
		H[0][k + 2] = v[2];

		x = H[k + 1][k];
		y = H[k + 2][k];

		if(k < n - 3)
		{
			z = H[k + 3][k];
		}
	}

	u[0] = x; 
	u[1] = y;

	b = house(u, v, 2, tol);

	// H[n - 1 : n, n - 2:n] = (I - b * v*v')*H[n - 1 : n, n - 2:n]
	// l[1:3] = v[1:2]' * H[n - 1 : n, n - 2 : n]
	for(i = 0; i < 3; i++)
	{
		l[i]  = v[0] * H[n - 2][n - 3 + i];
		l[i] += v[1] * H[n - 1][n - 3 + i];
	}
	// H[n - 1 : n, n - 2:n] = H[n - 1 : n, n - 2:n] - b * v[1:2] * l[1:3]
	for(i = 0; i < 3; i++)
	{
		H[n - 2][n - 3 + i] = H[n - 2][n - 3 + i] - b * v[0] * l[i];
		H[n - 1][n - 3 + i] = H[n - 1][n - 3 + i] - b * v[1] * l[i];
	}

	// Store esential part of u where new zeros where introduced.
	H[n - 2 + 1][n - 3] = v[1];
	H[n - 1 + 1][n - 3] = v[2];

	// H[1:n, n - 1:n] = H[1:n, n - 1:n] * (I - b * v*v')
	// l[1:n] = H[1:n, n - 1:n] * v[1:2]
	for(i = 0; i < n; i++)
	{
		l[i]  = H[i][n - 2] * v[0];
		l[i] += H[i][n - 1] * v[1]; 
	}
	// H[1:n, n - 1:n]  = H[1:n, n - 1:n]  - b * l[1:n] * v[1:2]'
	for(i = 0; i < n; i++)
	{
		H[i][n - 2] = H[i][n - 2] - b * l[i] * v[0];
		H[i][n - 1] = H[i][n - 1] - b * l[i] * v[1];
	}
	// Store esential part of u where new zeros where introduced.
	H[0][0 + k + 1] = v[0 + 1];
	H[0][1 + k + 1] = v[1 + 1];

	delete[] l;
}


void QR(Matrix& Q, Matrix& R)
{
	
}
