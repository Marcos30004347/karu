#include <cmath>
#include "algebra/matrix/Matrix.hpp"

namespace karu::algebra {

// Computes (a2 + b2)1/2 without destructive underflow or overflow.
float pythag(float a, float b)
{
	float absa, absb;

	absa = fabs(a);
	absb = fabs(b);

	if(absa > absb)
		return absa*sqrt(1.0+(absb/absa)*(absb/absa));

	return (absb == 0.0 ? 0.0 : absb * sqrt(1.0+(absa/absb)*(absa/absb)));
}

void svd(Matrix& M, Matrix& U, Matrix& D, Matrix& V) {

	int m = M.rows();
	int n = M.columns();

	U = Matrix(M);
	D = Matrix(n, 1);
	V = Matrix(n, n);

	int flag, i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z, *rv1, tmp;

	rv1 = new float[n]{0};
	g = scale = anorm = 0.0;

	// Householder reduction to bidiagonal form.
	for (i=1; i<=n; i++) 
	{
		l = i+1;
		rv1[i-1] = scale*g;
		g = s = scale = 0.0;

		if(i <= m) 
		{
			for (k=i; k<=m; k++)
			{
				scale += fabs(U[k-1][i-1]);
			}

			if(scale) 
			{
				for (k=i;k<=m;k++)
				{
					U[k-1][i-1] /= scale;
					s += U[k-1][i-1]*U[k-1][i-1];
				}
			
				f = U[i-1][i-1];
				g = -((f) >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				U[i-1][i-1] = f-g;
			
				for (j=l; j<=n; j++) 
				{
					s = 0.0;

					for (k=i; k<=m; k++)
					{
						s += U[k-1][i-1]*U[k-1][j-1];
					}
					
					f=s/h;
					
					for (k=i; k<=m; k++)
					{
						U[k-1][j-1] += f*U[k-1][i-1];
					}
				}
		
				for (k=i;k<=m;k++)
				{
					U[k-1][i-1] *= scale;
				}
			}
		}

		D[i-1][0] = scale*g;

		g = 0.0;
		s = 0.0;
		scale = 0.0;

		if (i <= m && i != n)
		{
			for (k=l; k<=n; k++)
			{
				scale += fabs(U[i-1][k-1]);
			} 
		
			if(scale)
			{
				for (k=l;k<=n;k++)
				{
					U[i-1][k-1] /= scale;
					s += U[i-1][k-1]*U[i-1][k-1];
				}

				f = U[i-1][l-1];
				g = -((f) >= 0.0 ? fabs(sqrt(s)) : -fabs(sqrt(s)));
				h = f*g-s;
				U[i-1][l-1] = f-g;
		
				for (k=l; k<=n; k++)
				{
				rv1[k-1]=U[i-1][k-1]/h;
				}
		
				for (j=l;j<=m;j++) 
				{
					for (s=0.0,k=l;k<=n;k++)
					{
						s += U[j-1][k-1]*U[i-1][k-1];
					}

					for (k=l;k<=n;k++)
					{
						U[j-1][k-1] += s*rv1[k-1];
					}
				}
				for (k=l;k<=n;k++)
				{
					U[i-1][k-1] *= scale;
				}
			}
		}

		tmp = (fabs(D[i-1][0])+fabs(rv1[i-1]));
		anorm = (anorm) > tmp ? (anorm) : tmp;
	}

	//Accumulation of right-hand transformations.
	for (i=n;i>=1;i--) 
	{
		if (i < n) 
		{
			if (g) 
			{
				//Double division to avoid possible underflow.
				for (j=l; j<=n; j++) 
				{
					V[j-1][i-1] = (U[i-1][j-1]/U[i-1][l-1])/g;
				}

				for (j=l;j<=n;j++) 
				{
					for (s=0.0,k=l;k<=n;k++)
					{
						s += U[i-1][k-1]*V[k-1][j-1];
					}
		
					for (k=l;k<=n;k++)
					{
						V[k-1][j-1] += s*V[k-1][i-1];
					}
				}
			}
			
			for (j=l;j<=n;j++) 
			{
				V[i-1][j-1]=V[j-1][i-1] = 0.0;
			}
		}
	
		V[i-1][i-1] = 1.0;
		g = rv1[i-1];
		l = i;
	}

	// Accumulation of left-hand transformations.
	for(i= (m < n ?m : n); i>=1; i--)
	{
		l = i+1;
		g = D[i-1][0];

		for (j=l;j<=n;j++)
		{
			U[i-1][j-1]=0.0;
		}

		if (g)
		{
			g=1.0/g;

			for (j=l;j<=n;j++) 
			{
				for (s=0.0,k=l;k<=m;k++)
				{
					s += U[k-1][i-1]*U[k-1][j-1];
				}

				f = (s/U[i-1][i-1])*g;

				for (k=i;k<=m;k++) 
				{
					U[k-1][j-1] += f*U[k-1][i-1];
				}
			}
		
			for (j=i;j<=m;j++)
			{
				U[j-1][i-1] *= g;
			} 
		} 
		else
		{
			for (j=i;j<=m;j++)
			{
				U[j-1][i-1]=0.0;
			}
		}
	
		U[i-1][i-1]++;
	}


	//Diagonalization of the bidiagonal form: Loop over
	for (k=n;k>=1;k--) 
	{ 
		//singular values, and over allowed iterations.
		for (its=1; its<=30; its++) 
		{
			flag=1;
				
			//Test for splitting.
			for (l=k;l>=1;l--)
			{ 
				nm = l-1; //Note that rv1[1] is always zero.
				if((float)(fabs(rv1[l-1])+anorm) == anorm)
				{
					flag=0;
					break;
				}
		
				if((float)(fabs(D[nm-1][0])+anorm) == anorm)
					break;
			}
			if(flag)
			{
				// Cancellation of rv1[l], if l > 1.
				c = 0.0;
				s = 1.0;
				for (i=l;i<=k;i++)
				{
					f = s*rv1[i-1];
					
					rv1[i-1] = c*rv1[i-1];
					
					if((float)(fabs(f)+anorm) == anorm) 
						break;
					
					g = D[i-1][0];
					
					h = pythag(f,g);
					
					D[i-1][0] = h;
					
					h = 1.0/h;
					c = g*h;
					s = -f*h;
					
					for(j=1;j<=m;j++)
					{
						y=U[j-1][nm-1];
						z=U[j-1][i-1];
						U[j-1][nm-1]=y*c+z*s;
						U[j-1][i-1]=z*c-y*s;
					}
				}
			}

			z = D[k-1][0];

			//Convergence.
			if (l == k) 
			{
				//Singular value is made nonnegative.
				if (z < 0.0)
				{
					D[k-1][0] = -z;
		
					for(j=1;j<=n;j++)
					{
						V[j-1][k-1] = -V[j-1][k-1];
					}
				}
				break;
			}

			if (its == 30)
			{
				fprintf(stderr, "SVD ERROR: no convergence in 30 svdcmp iterations\n");
				exit(1);
			} 
		
			//Shift from bottom 2-by-2 minor.
			x = D[l-1][0]; 
			nm = k-1;
			y = D[nm-1][0];
			g = rv1[nm-1];
			h = rv1[k-1];
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+((f) >= 0.0 ? fabs(g) : -fabs(g))))-h))/x;

			c = s = 1.0; //Next QR transformation:
			
			for(j=l; j<=nm; j++)
			{
				i = j+1;
				g = rv1[i-1];
				y = D[i-1][0];
				h = s*g;
				g = c*g;
				z = pythag(f,h);
				rv1[j-1] = z;
				c = f/z;
				s = h/z;
				f = x*c+g*s;
				g = g*c-x*s;
				h = y*s;
				y *= c;

				for (jj=1;jj<=n;jj++)
				{
					x = V[jj-1][j-1];
					z = V[jj-1][i-1];
					V[jj-1][j-1] = x*c+z*s;
					V[jj-1][i-1] = z*c-x*s;
				}

				z = pythag(f,h);
				D[j-1][0] = z; //Rotation can be arbitrary if z = 0.

				if (z) 
				{
					z=1.0/z;
					c=f*z;
					s=h*z;
				}

				f=c*g+s*y;
				x=c*y-s*g;

				for (jj=1;jj<=m;jj++)
				{
					y=U[jj-1][j-1];
					z=U[jj-1][i-1];
					U[jj-1][j-1]=y*c+z*s;
					U[jj-1][i-1]=z*c-y*s;
				}
			}

			rv1[l-1]=0.0;
			rv1[k-1]=f;
			D[k-1][0]=x;
		}
	}

	delete rv1;
}



}

