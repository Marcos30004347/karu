#pragma once

#include <cmath>
#include "algebra/matrix/Matrix.hpp"

namespace karu::algebra {

f32 norm2(f32* v, i32 n, f32 tol)
{
	f32 norm = 0;
	for(i32 i=0; i<n; i++)
		norm += v[i]*v[i];
	return sqrt(norm);
}


void golubKahanLanczosBidiagonalization(Matrix& A, f32* gamma, f32* phi, Matrix& u, Matrix& v)
{
	i32 i, j, k, m, n;
	
	m = A.rows();
	n = A.columns();

	v = Matrix(n,n);
	u = Matrix(m,n);

	phi[0] = 0;
	gamma[0] = 0;

	// choose v[0] so that it is unit 2-norm
	v[0][0] = 1;

	for(k=0; k<n; k++)
	{
		// u[k] = A*v[k] - b[k-1]*u[k-1]
		for(i=0; i<m; i++)
			for(j=0; j<n; j++)
				u[i][k] = u[i][k] + A[i][j]*v[j][k];
	
		if(k>0)
		for(i=0; i<m; i++)
				u[i][k] = u[i][k] - phi[k-1]*u[i][k-1];

		gamma[k] = 0;
		for(i=0; i<m; i++)
			gamma[k] = gamma[k] + (u[i][k] * u[i][k]);
		gamma[k] = sqrt(gamma[k]);

		for(i=0; i<m; i++)
			u[i][k] = u[i][k]/gamma[k];
		
		if(k < n-1)
		{
			for(i=0; i<n; i++)
				for(j=0; j<m; j++)
					v[i][k+1] = v[i][k+1] + A[j][i]*u[j][k];

			for(i=0; i<n; i++)
					v[i][k+1] = v[i][k+1] - gamma[k]*v[i][k];
		
			phi[k] = 0;
			for(i=0; i<n; i++)
				phi[k] = phi[k] + (v[i][k+1] * v[i][k+1]);
			phi[k] = sqrt(phi[k]);

			for(i=0; i<n; i++)
				v[i][k+1] = v[i][k+1]/phi[k];
		}
	}

	// for(i=0; i<n; i++)
	// 	std::cout << gamma[i] << " ";
	// std::cout << "\n";

	// for(i=0; i<n-1; i++)
	// 	std::cout << phi[i] << " ";
	// std::cout << "\n";

	// printMatrix(u);
	// std::cout << "\n";
	// printMatrix(v);

	// Matrix B(n,n, {
	// 	gamma[0], phi[0],    0,
	// 	0,    gamma[1], phi[1],
	// 	0,    0, 	    gamma[2],
	// });
	// std::cout << "\n";
	// std::cout << "\n";

	// printMatrix(u*B*transpose(v));
	// std::cout << "\n";
	// std::cout << "\n";

	// printMatrix(A - u*B*transpose(v));

}

}
