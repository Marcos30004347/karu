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

void golubKahanLanczosBidiagonalization(Matrix& A)
{
	f32 *b, *a;
	i32 i, j, k, m, n;
	
	m = A.rows();
	n = A.columns();

	Matrix v(n+1,n+1);
	Matrix u(m+1,n+1);

	a = new f32[n+1];
	b = new f32[n+1];

	b[0] = 0;
	a[0] = 0;

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
				u[i][k] = u[i][k] - b[k-1]*u[i][k-1];

		a[k] = 0;
		for(i=0; i<m; i++)
			a[k] = a[k] + (u[i][k] * u[i][k]);
		a[k] = sqrt(a[k]);

		for(i=0; i<m; i++)
			u[i][k] = u[i][k]/a[k];
		
		// v[i][k+1] = A'*u[i][k] - a[k]*v[k]
		for(i=0; i<n; i++)
			for(j=0; j<m; j++)
				v[i][k+1] = v[i][k+1] + A[j][i]*u[j][k];

		for(i=0; i<m; i++)
				v[i][k+1] = v[i][k+1] - a[k]*v[i][k];
	
		b[k] = 0;
		for(i=0; i<m; i++)
			b[k] = b[k] + (v[i][k+1] * v[i][k+1]);
		b[k] = sqrt(b[k]);

		for(i=0; i<m; i++)
			v[i][k+1] = v[i][k+1]/b[k];
	}

	for(i=0; i<n; i++)
		std::cout << a[i] << " ";
	std::cout << "\n";

	for(i=0; i<n-1; i++)
		std::cout << b[i] << " ";
	std::cout << "\n";

	printMatrix(u);
	printMatrix(v);

}

}
