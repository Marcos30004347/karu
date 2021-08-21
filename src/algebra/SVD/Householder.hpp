#pragma once

#include "algebra/matrix/Matrix.hpp"
#include "algebra/linear/Linear.hpp"

#include <cmath>

namespace karu::algebra {

f32 norm(f32* v, i32 n, f32 tol)
{
	f32 norm = 0;
	for(i32 i=0; i<n; i++)
		norm += v[i]*v[i];
	return sqrt(norm);
}


// f32 norm(f32* x, i32 n, f32 tol, f32 low, f32 high)
// {
// 	i32 i = 0;
// 	f32 sum = 0.0;
// 	f32 sc = 0.0;

// 	while(x[i++] <= tol) {}
// 	if(i >= n) return sum;

// 	while

// }

// Given a real n-vector x, the Householder matrix
// P = I - u*u'/tal can be constructed so that P is
// orthogonal and P*x = += |x| * [1, 0 ... 0]
// The function outputs z = Py given x
void householderOne(f32* x, f32* y, f32* z, i32 n)
{
	i32 i;
	f32 sxx, suy, s, tal, d; 

	f32* u = (f32*)malloc(sizeof(f32)*n);

	for(i=0; i<n; i++)
		sxx += x[i]*x[i];
	
	if(sxx == 0)
	{
		for(i=0; i<n; i++)
			z[i] = y[i];
	}
	else
	{
		s = sqrt(sxx);
		u[0] = x[0] + s;

		for(i=1; i<n; i++)
			u[i] = x[i];

		tal = sxx + s*x[0];

		for(i=0; i<n; i++)
			suy = suy + u[i]*y[i];

		d = suy/tal;
		
		for(i=0; i<n; i++)
			z[i] = y[i] - u[i]*d;
	}
	
	delete u;
}

// Given a real n-vector x, the Householder matrix
// P = I - u*u'/tal can be constructed so that P is
// orthogonal and P*x = += |x| * [1, 0 ... 0]
// The function outputs z = Py given x
void householderTwo(f32* x, f32* y, f32* z, i32 n)
{
	i32 i;
	f32 sxx, sxy, s, tal, d; 

	f32* u = (f32*)malloc(sizeof(f32)*n);

	for(i=0; i<n; i++)
		sxx += x[i]*x[i];
	
	if(sxx == 0)
	{
		for(i=0; i<n; i++)
			z[i] = y[i];
	}
	else
	{
		s = sqrt(sxx);

		u[0] = x[0] + s;
	
		for(i=0; i<n; i++)
			sxy = sxy + x[i]*y[i];
		
		z[0] = -sxy/s;
	
		d = (y[0] - z[0])/u[0];
	
		for(i=1; i<n; i++)
			z[i] = y[i] - x[i]*d;
	}
	
	delete u;
}

// Given x[0:n] this function computes a = alpha 
// and the vector u[0:n] such that (I - u*u')x = alpha [1,0...0]
void housegen(f32* x, f32* u, f32* a, i32 n, f32 tol)
{
	i32 i;

	f32 p, s, v;

	for(i=0; i<n; i++)
		u[i] = x[i];

	v = norm(x, n, tol);

	if(v == 0)
	{
		u[0] = sqrt(2);
		*a = v;
	} 
	else
	{
		if(u[0] >= 0.0)
		{
			p = u[0]/fabs(u[0]);
		}
		else
		{
			p = 1.0;
		}

		for(i=0; i<n; i++)
		{
			u[i] = (p/v)*u[i];
		}

		u[0] = 1 + u[0];

		s = sqrt(u[0]);
		
		for(i=0; i<n; i++)
		{
			u[i] = u[i]/s;
		}

		*a = -p * v;
	}
}


// Calculate X[:, k:n]*H given that H*z = -alhpa*e1
f32 applyHouseholderRight(Matrix& X, f32* z, i32 k, i32 n, i32 m, f32* Vk, f32 tol)
{
	i32 i, j, q;

	f32 alpha;

	f32* u = (f32*)malloc(sizeof(f32)*n - k);
	f32* b = (f32*)malloc(sizeof(f32)*m);

	housegen(z, u, &alpha, n - k, tol);

	// A*H = A(I - u*u') = A - (A*u)*u

	// We now have (I - u*u') = V

	// Vk = (I - u*u')'	
	for(i=0; i < n-k; i++)
	{
		for(j=0; j < n-k; j++)
		{
			if(i == j)
				Vk[i*(n-k) + j] = 1.;
			else
				Vk[i*(n-k) + j] = 0.;
		}
	}

	for(i=0; i < n-k; i++)
		for(j=0; j < n-k; j++)
			Vk[i*(n-k) + j] = Vk[i*(n-k) + j] - u[i]*u[j];
	

	// we need X[:][k:n]*V'
	
	// given that V is hermitian,
	// X[:][k:n]*V' 
	//	= X[:][k:n] * (I - u*u')
	//	= X[:][k:n] - (X[:][k:n]*u)*u'

	// X[:][k:n]*u
	for(i=0; i<m; i++)
	{
		b[i] = 0;
		for(j=k; j<n; j++)
		{
			b[i] += X[i][j]*u[j-k];
		}
	}

	//X[:][k:n] - b*u'
	for(i=0; i<m; i++)
	{
		for(j=k; j < n; j++)
		{
			X[i][j] = X[i][j] - b[i]*u[j-k];
			
			// Inverted results here
			// X[i][j] = -X[i][j];
		}
	}


	// Matrix zk(n-k, 1, z);
	// Matrix Vkk(n-k, n-k, Vk);
	// Vkk = Vkk*-1;
	// std::cout << "***\n";
	// std::cout << "alpha: "<< alpha << "\n\n";
	// std::cout << "P*z: [\n";
	// printMatrix(Vkk*zk);
	// std::cout << "]\n";
	delete u;

	return alpha;
}


void barlowBidiagonalization(Matrix X, f32* gamma, f32* phi, Matrix& u, Matrix& v, f32 tol, f32* p_)
{
	i32 i, j, k, q;
	i32 m = X.rows();
	i32 n = X.columns();
	
	f32 t, p;
	// f32* phi = (f32*)malloc(sizeof(f32)* n);
	// f32* gamma = (f32*)malloc(sizeof(f32)* n);

	f32* y = (f32*)malloc(sizeof(f32)* m);
	f32* z = (f32*)malloc(sizeof(f32)* m);
	f32* x = (f32*)malloc(sizeof(f32)* m);
	f32* Vk = (f32*)malloc(sizeof(f32)*m);

	Matrix tmp = Matrix(n, n);

	u = Matrix(m, n);
	v = Matrix(n, n);

	for(i=0;i<n;i++)
	{
		v[i][i] = 1.;
		tmp[i][i] = 1.;
	}
	
	for(i=0; i<m; i++)
	{
		y[i] = X[i][0];
	}
	
	gamma[0] = norm(y, m, tol);

	for(i=0; i<m; i++)
		u[i][0] = y[i]/gamma[0];

	q = n - 1;

	for(j=0; j<q; j++)
	{
		for(i=0; i<m; i++)
		{
			z[j] += X[i][j+1]*u[i][0];
		}
	}

	phi[1] = norm(z, q, tol);

	if(phi[1] != 0)
	{
		phi[1] = applyHouseholderRight(X, z, 1, n, m, Vk, tol);

		// accumulate V[k]
		// columns 0 is unnaltered
		for(j=0; j<1; j++)
		{
			for(i=0; i<n; i++)
			{
				tmp[i][j] = v[i][j]; 
			}
		}

		// V[:][1:n] = Vk*V[:][1:n]
		for(i=0; i<n; i++)
		{
			for(j=1; j<n; j++)
			{
				tmp[i][j] = 0.0;

				for(q=1; q<n; q++)
				{
					tmp[i][j] = tmp[i][j] + v[i][q]*Vk[(q - 1) * (n - 1) + j - 1];
				}
			}
		}

		v = tmp;
	}

	t = 0;

	for(k=1; k < n - 2; k++)
	{
	
		if(phi[k] == 0)
		{
			for(i=0;i<m; i++)
			{
				y[i] = X[i][k];
			}
		}
		else
		{
			for(i=0; i<m; i++)
			{
				y[i] = X[i][k] - (phi[k] * u[i][k-1]);
			}
		}

		gamma[k] = norm(y, m, tol);

		for(i=0; i < m; i++)
		{
			u[i][k] = y[i]/gamma[k];
		}

		q = n - k + 1;
	
		for(j=0; j<q; j++)
		{
			for(i=0; i<m; i++)
			{
				z[j] += X[i][k+1+j]*u[i][k];
			}
		}

		phi[k+1] = norm(z, q, tol);
	
		if(phi[k+1] != 0.0)
		{
			phi[k+1] = applyHouseholderRight(X, z, k + 1, n, m, Vk, tol);
			
			// Multiply V by {I(k), 0; 0, Vk

			// columns 0:k+1 is unnaltered
			for(j=0; j<k+1; j++)
			{
				for(i=0; i<n; i++)
				{
					tmp[i][j] = v[i][j]; 
				}
			}

			for(i=0; i<n; i++)
			{
				for(j=k+1; j<n; j++)
				{
					tmp[i][j] = 0.0;
					for(q=k+1; q<n; q++)
					{
						tmp[i][j] = tmp[i][j] + v[i][q]*Vk[(q - k + 1) * (n - 1) + j - k + 1];
					}
				}
			}

			v = tmp;
		}
		else 
		{
			phi[k+1] = 0;
		}

		t = fabs(phi[k+1]) + fabs(gamma[k+1]);
		if(t > p) p = t;
	}

	if(n > 2)
	{
		for(i=0; i < m; i++)
			phi[n-2] += u[i][n-3]*X[i][n-2]; 

		for(i=0; i < m; i++)
			y[i] = X[i][n-2] - phi[n-2]*u[i][n-3];
		
		gamma[n-2] = norm(y, m, tol);

		for(i=0;i<m;i++)
			u[i][n-2] = y[i]/gamma[n-2];
		
		t = fabs(phi[n-2]) + fabs(gamma[n-2]);
		if(t > p) p = t;
	}

	if(n > 1)
	{
		for(i=0;i<m;i++)
			phi[n-1] += u[i][n-2]*X[i][n-1]; 

		for(i=0;i<m;i++)
			y[i] = X[i][n-1] - phi[n-1]*u[i][n-2];
		
		gamma[n-1] = norm(y, m, tol);

		for(i=0;i<m;i++)
			u[i][n-1] = y[i]/gamma[n-1];

		t = fabs(phi[n-2]) + fabs(gamma[n-2]);
		if(t > p) p = t;
	}

	if(p_) *p_ = p;
	
	delete x;
	delete y;
	delete z;
	delete Vk;
}


// arrays e and q are gonna hold the values of the diagonals
// e is the upper diagonal and q the lower.
void golubReinschHouseholderBidiagonalization(
	Matrix& a,
	f32* q,
	f32* e,
	Matrix& u,
	Matrix& v,
	f32 tol = 2.22e-16,
	f32* x_ = nullptr
)
{
	i32 i, j, k, l, l1;
	f32 c, f, g, h, s, x, y, z;

	// m >= n assumed
	i32 m = a.rows();
	i32 n = a.columns();

	// Copy matrix a into u
	// u = Matrix(a);
	u = Matrix(m, m);

	for(i=0; i<m; i++)
		for(j=0; j<n; j++)
			u[i][j] = a[i][j];

	v = Matrix(n, n);

	g = 0;
	x = 0;

	// Householder reduction t o bidiagonal form
	for(i = 0; i < n; i++)
	{

		e[i] = g;

		s = 0;
		l = i+1;

		for(j = i; j < m; j++)
			s = s + (u[j][i] * u[j][i]);

		if(s < tol) 
			g = 0;
		else
		{
			f = u[i][i];
			g = f < 0 ? sqrt(s) : -sqrt(s);

			h = f * g - s;

			u[i][i] = f - g;

			for(j = l; j<n; j++)
			{
				s = 0;
			
				for(k = i; k<m; k++)
					s = s + (u[k][i] * u[k][j]);
			
				f = s / h;

				for(k=i; k<m; k++)
					u[k][j] = u[k][j] + (f * u[k][i]);
			}
		}
	
		q[i] = g; 
		s = 0;

		for(j = l; j < n; j++)
			s = s + (u[i][j] * u[i][j]);
		
		if(s < tol)
			g = 0;
		else
		{

			f = u[i][i+1];
			g = f < 0 ? sqrt(s) : -sqrt(s);
			
			h = f * g - s;
			u[i][i+1] = f - g;

			for(j = l; j<n; j++)
				e[j] = u[i][j] / h;
			
			for(j = l; j<m; j++)
			{
				s = 0;

				for(k = l; k < n; k++)
					s = s + (u[j][k] * u[i][k]);
				
				for(k = l; k < n; k++)
					u[j][k] = u[j][k] + (s * e[k]);
			}
		}
		y = fabs(q[i]) + fabs(e[i]);
		
		if(y > x) x = y;
	}

	// accumulation of right-hand transformations.
	for(i=n-1; i>=0; i--)
	{
		if(g != 0)
		{
			h = u[i][i+1] * g;
		
			for(j = l; j < n; j++)
				v[j][i] = u[i][j]/h;

			for(j = l; j < n; j++)
			{
				s = 0;
			
				for(k = l; k < n; k++)
					s = s + (u[i][k]*v[k][j]);
			
				for(k = l; k < n; k++)
					v[k][j] = v[k][j] + (s * v[k][i]);
			
			}
		}

		for(j=l; j<n; j++)
		{
			v[i][j] = 0;
			v[j][i] = 0;
		}

		v[i][i] = 1;

		g = e[i];

		l = i;
	}
	// accumulation of the left-hand transformations.
	for(i=n; i<m; i++)
	{
		for(j=n; j>m;)
			u[i][j] = 0;

		u[i][i] = 1;
	}

	for(i=n-1; i>=0; i--)
	{
		l = i + 1;
		g = q[i];

		for(j = l; j < m; j++)
			u[i][j] = 0;

		if(g != 0)
		{
			h = u[i][i] * g;
			for(j = l; j<m; j++)
			{
				s = 0;
				for(k = l; k<m; k++)
					s = s + (u[k][i] * u[k][j]);
				
				f = s / h;
				
				for(k = i; k < m; k++)
					u[k][j] = u[k][j] + (f * u[k][i]);
			}

			for(j=i; j<m; j++)
				u[j][i] = u[j][i]/g;
		}
		else 
		{
			for(j = i; j < m; j++)
				u[j][i] = 0;
		}
		
		u[i][i] = u[i][i] + 1;
	}

	if(x_) *x_ = x;
}

}
