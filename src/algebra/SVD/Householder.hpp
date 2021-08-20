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

// // Given x[0:n] this function computes a = alpha 
// // and the vector u[0:n] such that (I - u*u')x = alpha [1,0...0]
// void housegen(f32* x, f32* u, f32* a, i32 n, f32 tol)
// {
// 	i32 i;

// 	f32 p, v;

// 	for(i=0; i<n; i++)
// 		u[i] = x[i];
	

// 	v = norm(x, n, tol);

// 	if(v == 0)
// 	{
// 		u[0] = sqrt(2);
// 		*a = v;
// 		return;
// 	}

// 	if(u[0] >= 0.0)
// 	{
// 		p = u[0]/fabs(u[0]);
// 	}
// 	else
// 	{
// 		p = 1.0;
// 	}

// 	// WARNING: signal inversion wasent present
// 	// on the original algorithm and was put here
// 	// so that (I - u*u')x = alpha*e1 instead of
// 	// (I - u*u')x = -alpha*e1.
// 	for(i=0; i<n; i++)
// 		u[i] = (p/v)*u[i];

// 	u[0] = 1 + u[0];

// 	f32 s = sqrt(u[0]);

// 	for(i=0; i<n; i++)
// 		u[i] = u[i]/s;

// 	*a = -p * v;
// }


// Given x[0:n] this function computes a = alpha 
// and the vector u[0:n] such that (I - u*u')x = alpha [1,0...0]
void housegen(f32* x, f32* u, f32* a, i32 n, f32 tol)
{
	i32 i;

	f32 v;

	for(i=0; i<n; i++)
		u[i] = x[i];
	
	v = norm(x, n, tol);

	if(v == 0)
	{
		u[0] = sqrt(2);
		*a = v;
		return;
	}

	for(i=0; i<n; i++)
		u[i] = x[i]/v;

	if(u[0] >= 0.0)
	{
		u[0] = u[0] + 1.0;
		v = -v;
	}
	else
	{
		u[0] = u[0] - 1;
	}

	f32 s = sqrt(u[0]);

	for(i=0; i<n; i++)
		u[i] = u[i]/s;
	
	*a = v;
}

// Calculate X[:, k:n]*H given that H*z = -alhpa*e1
void applyHouseholderRight(Matrix& X, f32* z, i32 k, i32 n, i32 m, Matrix& Vk, f32 tol)
{
	i32 i, j, q;

	f32 alpha;

	f32* u = (f32*)malloc(sizeof(f32)*n - k);
	f32* b = (f32*)malloc(sizeof(f32)*m);

	housegen(z, u, &alpha, n - k, tol);

	// A*H = A(I - u*u') = A - (A*u)*u

	// We now have (I - u*u') = V

	// Vk = (I - u*u')'	
	Vk = Matrix(n-k, n-k);
	for(i=0; i < n-k; i++)
		Vk[i][i] = 1.;
	
	for(i=0; i < n-k; i++)
		for(j=0; j < n-k; j++)
			Vk[i][j] = Vk[i][j] - u[i]*u[j];
	

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
		for(j=k; j < n; j++)
			X[i][j] = X[i][j] - b[i]*u[j-k];


	Matrix zk(n-k, 1, z);
	std::cout << "***\n";
	printMatrix(Vk*zk);

	delete u;
}



void barlowBidiagonalization(Matrix& X, Matrix& u, Matrix& v, f32 tol = 2.22e-16)
{
	i32 i, j, k, q;
	i32 m = X.rows();
	i32 n = X.columns();
	
	f32* phi = (f32*)malloc(sizeof(f32)* n);
	f32* gamma = (f32*)malloc(sizeof(f32)* n);

	f32* y = (f32*)malloc(sizeof(f32)* m);
	f32* z = (f32*)malloc(sizeof(f32)* m);
	f32* x = (f32*)malloc(sizeof(f32)* m);

	u = Matrix(m, n);
	v = Matrix(n, n);

	Matrix Vk;

	for(i=0;i<n;i++)
		v[i][i] = 1.;
	
	for(i=0; i<m; i++)
		y[i] = X[i][0];
	
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
		applyHouseholderRight(X, z, 1, n, m, Vk, tol);
	}
	else
	{
		
		Vk = Matrix(n-1,n-1);
		
		for(i=0;i<n-1;i++)
			Vk[i][i] = 1.;
	}

	// return;

	for(k=1; k < n-1; k++)
	{
		std::cout << "phi: " << phi[k] << "\n";
	
		if(phi[k] == 0)
		{
			for(i=0;i<m; i++)
				y[i] = X[i][k];
		}
		else
		{
			for(i=0;i<m; i++)
				y[i] = X[i][k] - phi[k]*u[i][k-1];
		}

		gamma[k] = norm(y, m, tol);
		std::cout << "gamma[k]: " << gamma[k]  << "\n";

		for(i=0;i<m;i++)
			u[i][k] = y[i]/gamma[k];
		
		// for(i=0;i<m;i++)
		// 	std::cout << u[i][k] << " ";
		// std::cout << "\n";

		q = n - k + 1;
	
		for(j=0; j<q; j++)
		{
			for(i=0; i<m; i++)
			{
				z[j] += X[i][k+1+j]*u[i][0];
			}
		}
	
		// for(i=0;i<n - k + 1;i++)
		// 	std::cout << z[i] << " ";
		// std::cout << "\n";

		phi[k+1] = norm(z, q, tol);
	
		std::cout << "phi[k+1]: " << phi[k+1] << "\n";
		if(phi[k+1] != 0.0)
		{
			applyHouseholderRight(X, z, k + 1, n, m, Vk, tol);
		}
		else 
		{
			// V[k+1] = I(n-k)
			Vk = Matrix(n-k,n-k);

			for(i = 0;i < n-k; i++)
				Vk[i][i] = 1.;

			phi[k+1] = 0;
		}
	}

	if(n > 1)
	{
		for(i=0;i<m;i++)
			phi[n-1] += u[i][n-1]*X[i][n-1]; 

		for(i=0;i<m;i++)
			y[i] = X[i][n-1] - phi[n-1]*u[i][n-2];
		
		gamma[n-1] = norm(y, m, tol);

		for(i=0;i<m;i++)
			u[i][n-1] = y[n-1]/gamma[n-1];
	}

	delete x;
	delete y;
	delete z;
}

// arrays e and q are gonna hold the values of the diagonals
// e is the upper diagonal and q the lower.
void householderBidiagonalForm(
	Matrix& a,
	f32* q,
	f32* e,
	Matrix& u,
	Matrix& v,
	f32 tol = 2.22e-16,
	f32* c_ = nullptr, 
	f32* f_ = nullptr,
	f32* g_ = nullptr,
	f32* h_ = nullptr,
	i32* l_ = nullptr,
	f32* x_ = nullptr,
	f32* y_ = nullptr,
	f32* z_ = nullptr
)
{
	i32 i, j, k, l, l1;
	f32 c, f, g, h, s, x, y, z;

	// m >= n assumed
	i32 m = a.rows();
	i32 n = a.columns();

	// Copy matrix a into u
	u = Matrix(m, m);
	// u = Matrix(a);

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

	if(c_) *c_ = c;
	if(g_) *g_ = g;
	if(h_) *h_ = h;
	if(f_) *f_ = f;
	if(l_) *l_ = l;
	if(x_) *x_ = x;
	if(y_) *y_ = y;
	if(z_) *z_ = z;
}

}
