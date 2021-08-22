#pragma once 

#include <math.h>
#include "algebra/matrix/Matrix.hpp"

namespace karu::algebra  {


f32 sigmaMin(f32* s, f32* e, i32 m, i32 n)
{
	i32 j;
	f32 *l, *u, b_inf, b_one, sigma;
	
	l = (f32*)malloc(sizeof(f32)*n);
	u = (f32*)malloc(sizeof(f32)*n);

	l[n-1] = fabs(s[n-1]);

	for(j = n-2; j >= 0; j--)
		l[j] = fabs(s[j]) * (l[j+1] / (l[j+1] + fabs(e[j])));

	u[0] = fabs(s[0]);

	for(j=0; j<n-1; j++)
		u[j+1] = fabs(s[j+1]) * (u[j] / (u[j] + fabs(e[j])));

	b_inf = l[0];
	b_one = u[0];

	for(j=1; j<n; j++)
		b_inf = std::min(b_inf, l[j]);

	for(j=1; j<n; j++)
		b_one = std::min(b_one, u[j]);

	sigma = std::min(b_one, b_inf);

	delete l;
	delete u;

	return sigma;
}


void bidiagonalSVD(f32* s, f32* e, i32 m, i32 n, i32 maxit, f32 tol)
{
	i32 i, j, i_, i_min, i_max;
	f32 sig_min, sig_max, thresh;

	sig_min = sigmaMin(s, e, m, n);

	thresh = tol * sig_min;

	for(i=0; i<n; i++)
		if(e[i] < thresh) e[i] = 0;

	while(true)
	{
		// Find bottommost nonscalar unreduced block
		// diagonal submatrix of B
		i_max = -1;
		for(i=0; i<n-1; i++)
			if(fabs(e[i]) > thresh)
			{
				i_max = i;
				break;
			}		
	
		if(i_max == -1)
			i_max = n-1;

		if(i_max == 1)
		{
			// TODO: goto done
		}

		i_ = -1;
		for(i=n-2; i>=0; i--)
			if(fabs(e[i]) <= thresh)
			{
				i_ = i;
				break;
			}		
		if(i_ == -1)
			i_ = 0;

		i_min = i_ + 1;

		// Apply algorithm to unreduced block diagonal 
		// submatrix from i_min to i_max
		if(i_max == i_min + 1)
		{
			// 2x2 submatrix handle specially
			// compute 2x2 svd submatrix, setting ei to 0
			continue;
		}
	}
}

}
