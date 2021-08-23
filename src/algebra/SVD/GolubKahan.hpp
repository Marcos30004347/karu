#pragma once

#include "algebra/matrix/Matrix.hpp"
#include <cmath>
#include <iomanip>
namespace karu::algebra {

i32 sign(f32 v)
{
	if(v > 0) return 1;
	if(v == 0) return 0;
	return -1;
}

// GvL pg. 216 : algo 5.1.3
void givens(f32 a, f32 b, f32* c, f32* s, f32 tol)
{
	// Computes scalars c and s such that
	//   [c, s; -s, c].T * [a, b] = [r, 0]
	f32 r;

	if(fabs(b) <= tol)
	{
		*c = 1;
		*s = 0;
	}
	else
	{
		// r = hypot(y, z); // safe for underflow and overflow
		
		// c = y/r;
		// s = z/r;
		if(fabs(b) > fabs(a))
		{
			r = -a / b;
			*s = 1 / sqrt(1 + (r * r));
			*c = *s * r;
		}
		else
		{
			r = -b / a;
			*c = 1 / sqrt(1 + (r * r));
			*s = *c * r;
		}
	}
}

// GvL pg 216 section 5.1.9
// Computes [c, s; -s,c].T * A[:, [i, k]]
void leftGivens(Matrix& A, f32 c, f32 s, i32 i, i32 k, f32 tol)
{
	// Apply [c, s; -s, c].T to A[[i,k], :]
	i32 j, n, m;
	f32 t1, t2;

	m = A.rows();
	n = A.columns();

	for(j=0; j<n; j++)
	{
		t1 = A[i][j];
		t2 = A[k][j];

		A[i][j] = c * t1 - s * t2;
		A[k][j] = s * t1 + c * t2;

		if(fabs(A[i][j]) <= tol)
			A[i][j] = 0.0;

		if(fabs(A[k][j]) <= tol)
			A[k][j] = 0.0;
	}
}



// Modified leftGivens(GvL pg 216 section 5.1.9) to work on
// diagonal matrice only, d is an array containing 5 diagonals
// of the matrice where the main diagonal is d[2], d[1]
// is the diagonal above d[2] and d[0] is above d[1], 
// d[3] is the super diagonal above d[2] and d[4] is the
// diagonal above d[3], d[2][0] is the first element of the matrix,
// elements of the diagonals that should be above it should equal zero
// and elements that come after d[2][n-1](last element of the matrix)
// should also be zero
void diagLeftGivens(f32* d[5], i32 m, i32 n, f32 c, f32 s, i32 i, i32 k, f32 tol)
{

	i32 col = 0;
	f32 t1, t2;

	i = (i + 2); // add main diagonal start idx
	k = (k + 2); // add main diagonal start idx

	while(i > 4)
	{
		i = i - 1;
		k = k - 1;
		col = col + 1;
	}

	if(i == 4 || k > 4)
	{
		t1 = d[i][col];
		t2 = 0;
	
		d[i][col] = c * t1 - s * t2;

		if(fabs(d[i][col]) <= tol)
			d[i][col] = 0.0;

		i = i - 1;
		k = k - 1;
	
		col++;
	}
	// d[0][0],       0,       0,
	// d[1][0], d[0][1], 			 0,
	// d[2][0], d[1][1], d[0][2],
	// d[3][0], d[2][1], d[1][2],
	// d[4][0], d[3][1], d[2][2]
	//      0,  d[4][1], d[3][2]
	// 			0,				0, d[4][2]

	for(; i >= 0; col++)
	{
		t1 = d[i][col];
		t2 = d[k][col];
		
		d[i][col] = c * t1 - s * t2;
		d[k][col] = s * t1 + c * t2;

		if(fabs(d[i][col]) <= tol)
			d[i][col] = 0.0;

		if(fabs(d[k][col]) <= tol)
			d[k][col] = 0.0;

		i = i - 1;
		k = k - 1;
	}
}


// GvL pg 216 section 5.1.9
// Computes A[:, [i, k]] * [c, s; -s,c]
void rightGivens(Matrix& A, f32 c, f32 s, i32 i, i32 k, f32 tol)
{
	// Apply [c, s; -s, c].T to A[[i,k], :]
	i32 j, n, m;
	f32 t1, t2;

	m = A.rows();
	n = A.columns();

	for(j=0; j<m; j++)
	{
		t1 = A[j][i];
		t2 = A[j][k];

		A[j][i] = c * t1 - s * t2;
		A[j][k] = s * t1 + c * t2;

		if(fabs(A[j][i]) <= tol)
			A[j][i] = 0.0;

		if(fabs(A[j][k]) <= tol)
			A[j][k] = 0.0;

	}
}


// Modified rightGivens(GvL pg 216 section 5.1.9) to work on
// diagonal matrice only, d is an array containing 5 diagonals
// of the matrice where the main diagonal is d[2], d[1]
// is the diagonal above d[2] and d[0] is above d[1], 
// d[3] is the super diagonal above d[2] and d[4] is the
// diagonal above d[3], d[2][0] is the first element of the matrix,
// elements of the diagonals that should be above it should equal zero
// and elements that come after d[2][n-1](last element of the matrix)
// should also be zero
void diagRightGivens(f32* d[5], i32 m, i32 n, f32 c, f32 s, i32 i, i32 k, f32 tol)
{
	i32 j, t, q;
	f32 t1, t2;

	t1 = d[0][i];
	t2 = 0;
	
	d[0][i] = c * t1 - s * t2;

	for(q=0; q<4; q++)
	{
		t1 = d[q+1][i];
		t2 = d[q][k];

		d[q+1][i] = c * t1 - s * t2;
		d[q][k] = s * t1 + c * t2;

		if(fabs(d[q+1][i]) <= tol)
		{
			d[q+1][i] = 0.0;
		}

		if(fabs(d[q][k]) <= tol)
		{
			d[q][k] = 0.0;
		}
	}

	t1 = 0;
	t2 = d[4][k];

	d[4][k] = s * t1 + c * t2;

	if(fabs(d[4][k]) <= tol)
	{
		d[4][k] = 0.0;
	}
}

// alocate 3 extra diagonals for B and set the values to zero
// the B = diags(diag[0], diag[1], diag = diag[2], sdiag = diag[3], diag[3], diag[4])
f32** buildDiagMatrix(f32* diag, f32* sdiag, i32 m , i32 n)
{
	f32** diags = new f32*[5];

	diags[0] = new f32[n];
	diags[1] = sdiag;
	diags[2] = diag;
	diags[3] = new f32[n];
	diags[4] = new f32[n];

	std::fill(&diags[0][0], &diags[0][n-1] + 1, 0);
	std::fill(&diags[3][0],  &diags[3][n-1]  + 1, 0);
	std::fill(&diags[4][0], &diags[4][n-1] + 1, 0);
	
	return diags;
}

// delete extra diagonals in diags
// pointer diags[2] and diags[3]
// stays allocated
void freeDiagMatrix(f32** diags)
{
	delete diags[0];
	delete diags[3];
	delete diags[4];
	delete[] diags;
}

f32 trailing2x2Eigenvalue(f32* b_diag, f32* b_sdiag, i32 m, i32 n, f32 tol)
{
	f32 a, b, c, d, mu, t11, t12, t21, t22;
	// 1. Find y and z 
	//	 Let u be the eigenvalue of the trailing 2x2 submatrix
	// 	 of T = B'*B that is closer to t_nn
	
	//	T[0:2, 0:2] = B[:,-2:].T*B[:,-2:]
	//  B[:, -2] = {{0,..., 0, a, c}, { 0,..., 0, b, d}}.T
	a = b_sdiag[n - 2];
	b = b_sdiag[n - 1];
	c = b_diag[n - 2];
	d = b_diag[n - 1];

	t11 = a*a + c*c;
	t12 = c*b;
	t21 = b*c;
	t22 = b*b + d*d;

	// std::cout << a << " " << b << "\n";
	// std::cout << c << " " << d << "\n";
	// std::cout << "\n";
	// std::cout << t11 << " " << t12 << "\n";
	// std::cout << t21 << " " << t22 << "\n";

	// compute wilkinson shift value "µ" ("mu")
	d = (t11 - t22) / 2.0;
	mu = t22 - (sign(d) * (t21 * t21)) / (fabs(d) + sqrt((d * d) + (t21 * t21)));
	// mu = t22 - (t21 * t21) / (d * sign(d) * sqrt((d * d) + (t21 * t21)));
	
	return mu;
}

void debugPrintDiagMatrix(f32** diag, i32 m, i32 n)
{
	// d[0][0],       0,       0,
	// d[1][0], d[0][1], 			 0,
	// d[2][0], d[1][1], d[0][2],
	// d[3][0], d[2][1], d[1][2],
	// d[4][0], d[3][1], d[2][2]
	//      0,  d[4][1], d[3][2]
	// 			0,				0, d[4][2]
	
	i32 i,j;
	for(i = 0; i < n; i++)
	{
		i32 d = i + 2;

		for(j = d; d - j < n && j >= 0; j--)
		{
			if(j > 4)
			{
				std::cout << std::fixed << std::setprecision(5) << std::setw(8) << 0.000 << ", ";
				continue;
			}

			std::cout << std::fixed << std::setprecision(5) << std::setw(8) << diag[j][d - j] << ", ";
		}
	
		d = d - j;
	
		while(d <= 3)
		{
			std::cout << std::fixed << std::setprecision(5) << std::setw(8) << 0.000 << ", ";
			d++;
		}
		std::cout << "\n";
	}

	for(i = 0; i < m - n; i++)
	{
		for(j=0; j < n; j++)
		{
			std::cout << std::fixed << std::setprecision(5) << std::setw(8) << 0.000 << ", ";
		}
		std::cout << "\n";
	}
}

void golubKahanStep(f32* b_diag, f32* b_sdiag, i32 m, i32 n, Matrix &uT, Matrix &v, f32 tol)
{
	// B[0:m-1, 0:n-1] = bidiag(b_diag, b_sdiag)
	i32 i, k;
	f32 y, z, c, s, mu, t11, t12;
	
	// uT = identity(m,m); 
	// v  = identity(n,n); 

	t11 = b_diag[0] * b_diag[0];
	t12 = b_diag[0] * b_sdiag[1];

	// std::cout << t11 << " " << t12  << "\n";

	mu = trailing2x2Eigenvalue(b_diag, b_sdiag, m, n, tol);

	// std::cout << "mu: " << mu << "\n";

	y = t11 - mu;
	z = t12;
	
	f32** diags = buildDiagMatrix(b_diag, b_sdiag, m, n);

		// Debug
		// std::cout << "\n";
		// debugPrintDiagMatrix(diags, m, n);
		// std::cout << "\n";
		
	bool gather_rotations = true;

	// for(k = from /* 0 */ ; k < to /* n - 1 */; k++)
	for(k = 0; k < n - 1; k++)
	{
		givens(y, z, &c, &s, tol);
		diagRightGivens(diags, m, n, c, s, k, k + 1, tol);
	
			// Debug
			// std::cout << "B * G(" << k << ", " << k+1 << ", theta)\n";
			// debugPrintDiagMatrix(diags, m, n);

		y = diags[2][k]; z = diags[3][k];

		// if(gather_rotations)
		// {
		rightGivens(v, c, s, k, k + 1, tol);
		// }

		// gather_rotations = true;
	
		givens(y, z, &c, &s, tol);
		diagLeftGivens(diags, m, n, c, s, k, k + 1, tol);

			// Debug
			// std::cout << "G(" << k << ", " << k+1 << ", theta).T * B\n";
			// debugPrintDiagMatrix(diags, m, n);

		// if(gather_rotations)
		// {
		// Save U.T
		leftGivens(uT, c, s, k, k+1, tol);
		// }

		if(k < n - 2)
		{
			y = diags[1][k+1]; z = diags[0][k+2];
		}
	}
	// std::cout << std::scientific << std::numeric_limits<f32>::epsilon() << "\n";
	// std::cout << "\n";	
	// std::cout << "uT:\n";
	// printMatrix(uT);
	// std::cout << "\n";
	// std::cout << "v:\n";
	// printMatrix(v);
	// std::cout << "\n";
	
	// Matrix b(4,4, {
	// 	1, 1, 0, 0,
	// 	0, 2, 1, 0,
	// 	0, 0, 3, 1,
	// 	0, 0, 0, 4.
	// });
	// std::cout << "b:\n";
	// printMatrix(b);
	// std::cout << "\n";
	// printMatrix(transpose(uT)*uT*b*v*transpose(v));
	// std::cout << "\n";
	// debugPrintDiagMatrix(diags, m, n);


	// for(i=0; i<n; i++)
	// 	std::cout << b_sdiag[i] << " ";
	// std::cout << "\n";
	// for(i=0; i<n; i++)
	// 	std::cout << b_diag[i] << " ";
	// std::cout << "\n";


	// freeDiagMatrix(diags);

}

void findLimits(f32* b_diag, f32* b_sdiag, i32 m, i32 n, i32* p, i32 *q, f32 tol)
{
	i32 i;
	
	*q = -1;
	*p = 0;

	for(i = n - 1; i>=1; i--)
	{
		// pick largest q such that B33 is diagonal
		if(*q == -1 && fabs(b_sdiag[i]) >= tol)
		{
			*q = n - (i + 1);
		}

		// pick smalest p such that B22 has nonzero on superdiagonal
		if(*q != -1 && fabs(b_sdiag[i]) <= tol)
		{
			*p = i;
			break;
		}
	}

	// bool found_start = false;
	// i32 i;

	// for(i = n-1; i>=0; i--)
	// {
	// 	if(found_start)
	// 	{
	// 		if(fabs(b_sdiag[i]) <= tol)
	// 		{
	// 			*p = i;
	// 			break;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		if(fabs(b_sdiag[i]) <= tol)
	// 		{
	// 			*p = *q = i;
	// 			found_start = true;
	// 		}
	// 	}
	// }

	// if(!found_start)
	// {
	// 	*p = 0;
	// 	*q = 0;
	// }

	// *q = *q + 1;
}

void clearDiag(f32* b_diag, f32* b_sdiag, i32 m, i32 n, f32 tol)
{
	i32 i;
	for(i = 0; i < n; i++)
	{
		if(fabs(b_diag[i]) <= tol)
			b_diag[i] = 0.0;
		if(fabs(b_sdiag[i]) <= tol)
			b_sdiag[i] = 0.0;
	}
}

void golubKahanSVD(f32* b_diag, f32* b_sdiag, i32 m, i32 n, Matrix &uT, Matrix &v, f32 tol)
{
	i32 i, j, q = 0, p = -1;
	
	while(q != n - q)
	{
		std::cout << p << " " << n - p - q << " " << q << "\n";

		// Remove elements that are less that the tolerance
		clearDiag(b_diag, b_sdiag, m, n, tol);
		
		// Find the largest q and the smallest p such that if
		// B = [B11, 0, 0; 0, B22, 0; 0, 0, B33] then B33 is diagonal
		// and B22 has nonzero superdiagonal
		findLimits(b_diag, b_sdiag, m, n, &p, &q, tol);

		// B11 = B[0:p-1, 0:p-1]
		// B22 = B[p:n - q - 1, p:n - q - 1]
		// B33 = B[n - q:n, n - q:n]

		std::cout << p << " " << n - p - q << " " << q << "\n";

		f32** b = buildDiagMatrix(b_diag, b_sdiag, m, n);
		debugPrintDiagMatrix(b, m, n);
		freeDiagMatrix(b);

		if(q < n)
		{
			// if any entry in B22 is zero, then zero
			// the superdiagonal entry in the same row.
			bool zeroed = false;

			for(i=p; i < n - p - (q - 1); i++)
			{
				if((i + 1 < n) && fabs(b_diag[i]) <= tol)
				{
					zeroed = true;
					b_sdiag[i+1] = 0.0;
				}
			}

			if(!zeroed)
			{
				f32* b22_diag  = new f32[n - p - q];
				f32* b22_sdiag = new f32[n - p - q];

				// Matrix ut_ = identity(n - p - q, n - p - q);
				// Matrix v_  = identity(n - p - q, n - p - q);

				for(i = p; i < n - p - q; i++)
				{
					b22_diag[i - p] = b_diag[i];
					b22_sdiag[i - p] = b_sdiag[i];
				}
		
				golubKahanStep(b22_diag, b22_sdiag, n - p - q, n - p - q, uT, v, tol);

				for(i = p; i < n - p - q; i++)
				{
					b_diag[i] = b22_diag[i - p];
					b_sdiag[i] = b22_sdiag[i - p];
				}

				// delete b22_diag;
				// delete b22_sdiag;
		
			}
		}
	}

}

}
