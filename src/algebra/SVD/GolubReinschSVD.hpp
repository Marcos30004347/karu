#pragma once

#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"
#include "algebra/SVD/Householder.hpp"

#include <cmath>

namespace karu::algebra {

f32 pythag(f32 a,f32 b)
{
	f32 c = fabs(a);
	f32 d = fabs(b);
	if (c > d)
		return c*sqrt(1.0 + ((d/c) * (d/c)));
	else
	{
			if(d == 0.0) return 0.0;
			else return d*sqrt(1.0 + ((c/d) * (c/d)));
	}
}
// Handbook Series Linear Algebra
// Singular Value Decomposition and Least Squares Solutions*
// Contributed by G. H. GOLUB and C. REINSCH
void golubReinschBidiabonalSVD(f32* gamma, f32*phi, i32 m, i32 n, Matrix& u, f32* q, Matrix& v, f32 tol, f32 eps, f32 x)
{
	i32 it, i, j, k, l, l1;

	f32 s, c, f, g, h, y, z;

	bool goto_test_f_convergence;
	
	for(i=0; i<n; i++)
	{
		q[i] = gamma[i];
	}

	eps = eps * x;
	for(k = n-1; k >= 0; k--)
	{
		while(true)
		{
			// test f splitting:
			goto_test_f_convergence = false;
			for(l=k; l>=0; l--)
			{
				if(fabs(phi[l]) <= eps)
				{
					goto_test_f_convergence = true;
					break;
				}
				if(fabs(q[l-1]) <= eps)
				{
					break;
				}
			}
		
			// cancelation
			if(!goto_test_f_convergence)
			{
				// cancelation of e[l] if l > 0
				c = 0.0;
				s = 1.0;

				l1 = l - 1;

				for(i=l; i<=k; i++)
				{
					f = s * phi[i];

					phi[i] = c * phi[i];

					if(fabs(f) <= eps) break;
	
					g = q[i];

					h = pythag(f,g);

					q[i] = h;

					c = g/h;
					s = -f/h;

					for(j=0; j<m; j++)
					{
						y = u[j][l1];
						z = u[j][i];
						u[j][l1] = y * c + z * s;
						u[j][i] = -y * s + z * c;
					}
				}
			}

			// test f convergence:
			z = q[k];

			if(l == k)
			{
				// convergence
				if(z < 0.0)
				{
					q[k] = -z;
					for(j=0; j<n; j++)
						v[j][k] = -v[j][k];
				}
	
				break;
			}

      // shift from bottom 2x2 minor
			x = q[l];
			y = q[k-1];
			g = phi[k-1];
			h = phi[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);

			if(f < 0)
			{
				f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x;
			}
			else 
			{
				f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x;
			}
	
			c = 1.0;
			s = 1.0;

			for(i=l+1; i<=k; i++)
			{
				g = phi[i];
				y = q[i];

				h = s * g;
				g = c * g;

				z = pythag(f, h);

				phi[i-1] = z;

				c = f / z;
				s = h / z;

				f = x * c + g * s;
				g = -x * s + g * c;
				h = y * s;
				y = y * c;
				

				for(j=0; j < n; j++)
				{
					x = v[j][i-1];
					z = v[j][i];
					v[j][i-1] = x * c + z * s;
					v[j][i] = -x * s + z * c;
				}
			
				z = pythag(f, h);
				q[i-1] = z;
				c = f / z;
				s = h / z;
				f = c * g + s * y;
				x = -s * g + c * y;
		
				for(j = 0; j<m; j++)
				{
					y = u[j][i-1];
					z = u[j][i];

					u[j][i-1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				}
			}

			phi[l] = 0;
			phi[k] = f;
			q[k] = x;
		}
	}
}

}
