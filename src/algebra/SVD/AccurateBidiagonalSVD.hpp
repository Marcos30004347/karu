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

}
