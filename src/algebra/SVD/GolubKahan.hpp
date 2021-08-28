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

void sortUV(f32* s, u32 n, Matrix& u, Matrix& v)
{
	i32 i, k, rows, cols, j, i_last;
	f32 s_last, tmp;

	for (i = 0; i < n; ++i) 
	{
		s_last = s[i];
		i_last = i;
	
		for (j = i + 1; j < n; ++j)
		{
			if (fabs(s[j]) >  fabs(s_last))
			{
				s_last = s[j];
				i_last = j;
			}
		}
		
		if (i_last != i)
		{
			tmp;
			tmp = s[i];

			s[i] = s[i_last];
			s[i_last] = tmp;

			rows = u.rows();
			cols = u.columns();
					
			for (k = 0; k < rows; ++k) 
			{
				tmp = u[k][i];
				u[k][i] = u[k][i_last];
				u[k][i_last] = tmp;
			}

			rows = v.rows();
			cols = v.columns();
		
			for (k = 0; k < rows; ++k) 
			{
				tmp = v[i][k];
				v[i][k] = v[i_last][k];
				v[i_last][k] = tmp;
			}


		}
	}
}

// GvL pg. 216 : algo 5.1.3
void givens(f32 a, f32 b, f32* c, f32* s, f32 tol)
{
	// Computes scalars c and s such that
	//   [c, s; -s, c].T * [a, b] = [r, 0]
	// if(fabs(a) <= tol) a = 0.0;
	
	if(b == 0)
	{
		*c = 1.0;
		*s = 0.0;
	}
	else
	{
		f32 r = hypot(a, b);
		*c = a/r;
		*s = -b/r;
	}

	// return;

	// if(fabs(b) <= tol)
	// {
	// 	*c = 1;
	// 	*s = 0;
	// }
	// else
	// {
	// 	// r = hypot(y, z); // safe for underflow and overflow
		
	// 	// c = y/r;
	// 	// s = z/r;
	// 	if(fabs(b) > fabs(a))
	// 	{
	// 		r = -a / b;
	// 		*s = 1 / sqrt(1 + (r * r));
	// 		*c = *s * r;
	// 	}
	// 	else
	// 	{
	// 		r = -b / a;
	// 		*c = 1 / sqrt(1 + (r * r));
	// 		*s = *c * r;
	// 	}
	// }
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

		// if(fabs(A[i][j]) <= tol)
		// 	A[i][j] = 0.0;

		// if(fabs(A[k][j]) <= tol)
		// 	A[k][j] = 0.0;
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
	// std::cout << i << " " << k << "\n";

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

		// std::cout << i << " " << "*" << " " << col << "\n";
		
		i = i - 1;
		k = k - 1;
	
		col++;
	}
	// d[0][0],       0,       0, 			0,
	// d[1][0], d[0][1], 			 0, 			0,


	// d[2][0], d[1][1], d[0][2], 			0,
	// d[3][0], d[2][1], d[1][2], d[0][3],
	// d[4][0], d[3][1], d[2][2], d[1][3],
	//      0,  d[4][1], d[3][2], d[2][3],


	// 			0,				0, d[4][2],	d[3][3],

	for(; col < n && i >= 0; col++)
	{
		// std::cout << i << " " << k << " " << col << "\n";
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

	if(col < n)
	{
		t1 = 0;
		t2 = d[k][col];
		
		d[k][col] = s * t1 + c * t2;

		if(fabs(d[k][col]) <= tol)
			d[k][col] = 0.0;
		
		// std::cout << "*" << " " << k << " " << col << "\n";
	}

	// std::cout << "\n";
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

		// if(fabs(A[j][i]) <= tol)
		// 	A[j][i] = 0.0;

		// if(fabs(A[j][k]) <= tol)
		// 	A[j][k] = 0.0;
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

	for(q=0; q < 4; q++)
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

	std::fill(diags[0], diags[0] + n, 0);
	std::fill(diags[3], diags[3] + n, 0);
	std::fill(diags[4], diags[4] + n, 0);
	
	return diags;
}

// delete extra diagonals in diags
// pointer diags[2] and diags[3]
// stays allocated
void freeDiagMatrix(f32** diags)
{
	diags[1] = nullptr;
	diags[2] = nullptr;
	
	delete[] diags[0];
	delete[] diags[3];
	delete[] diags[4];
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

	// compute wilkinson shift value "µ" ("mu")
	d = (t11 - t22) / 2.0;
	mu = t22 - (sign(d) * (t21 * t21)) / (fabs(d) + sqrt((d * d) + (t21 * t21)));
	
	return mu;
}


f32 trailing2x2Eigenvalue(Matrix& B, i32 p, i32 q, f32 tol)
{
	i32 m,n;
	f32 a, b, c, d, mu, t11, t12, t21, t22;
	// 1. Find y and z 
	//	 Let u be the eigenvalue of the trailing 2x2 submatrix
	// 	 of T = B'*B that is closer to t_nn
	
	//	T[0:2, 0:2] = B[:,-2:].T*B[:,-2:]
	//  B[:, -2] = {{0,..., 0, a, c}, { 0,..., 0, b, d}}.T
	// std::cout << p << " " << q << "\n";

	// printSubMatrix(B, p, q, p, q);

	if(q >= 3)
	{
		a = B[q - 3][q - 2];
		b = B[q - 2][q - 1];
		c = B[q - 2][q - 2];
		d = B[q - 1][q - 1];
	}
	else 
	{
		a = B[q - 2][q - 2];
		b = B[q - 2][q - 1];
		c = B[q - 1][q - 2];
		d = B[q - 1][q - 1];
	}

	// std::cout << a << " " << b << " " << c << " " << d << "\n";

	t11 = a*a + c*c;
	t12 = c*b;
	t21 = b*c;
	t22 = b*b + d*d;
	// std::cout << t11 << " " << t12 << " " << t21 << " " << t22 << "\n";

	// compute wilkinson shift value "µ" ("mu")
	d = (t11 - t22) / 2.0;
	mu = t22 - (sign(d) * (t21 * t21)) / (fabs(d) + sqrt((d * d) + (t21 * t21)));
	// std::cout << mu << "\n";
	
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
				std::cout << std::scientific << std::setprecision(5) << std::setw(12) << 0.000 << ", ";
				continue;
			}

			std::cout << std::scientific << std::setprecision(5) << std::setw(12) << diag[j][d - j] << ", ";
		}
	
		d = d - j;
	
		while(d < n)
		{
			std::cout << std::scientific << std::setprecision(5) << std::setw(12) << 0.000 << ", ";
			d++;
		}
		std::cout << "\n";
	}

	for(i = 0; i < m - n; i++)
	{
		for(j=0; j < n; j++)
		{
			std::cout << std::scientific << std::setprecision(5) << std::setw(12) << 0.000 << ", ";
		}
		std::cout << "\n";
	}
}

void golubKahanStep(Matrix& B, i32 p, i32 q, Matrix &uT, Matrix &v, f32 tol)
{
	i32 i, k;
	f32 y, z, c, s, mu, t11, t12;

	t11 = B[p][p] * B[p][p];
	t12 = B[p][p] * B[p][p+1];
	mu = trailing2x2Eigenvalue(B, p, q, tol);
	y = t11 - mu;
	z = t12;

	// std::cout << "mu: " << mu << "\n";

	// std::cout << "###########S##########" << "\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "#####################" << "\n";

	givens(y, z, &c, &s, tol);
	rightGivens(B, c, s, p, p + 1, tol);
	leftGivens(v, c, s, p, p + 1, tol);

	// std::cout << "ASASSA" << "\n";

	// std::cout << "#####################" << "\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "#####################" << "\n";

	for(k = p; k < q - 2; k++)
	{
		// std::cout << "###########HIHI##########" << "\n";
	
		// std::cout << c << " " << s << " " << k << " " << k+1 << "\n";
		// std::cout << "V\n";
		// printMatrix(v);
		// std::cout << "IS\n";
		// printMatrix(uT*B*v);
		// std::cout << "B Matrix:\n";
		// printSubMatrix(B, p, q, p, q);
		// std::cout << "\n";
	
		y = B[k][k]; z = B[k+1][k];

		givens(y, z, &c, &s, tol);
		leftGivens(B, c, s, k, k + 1, tol);
		rightGivens(uT, c, s, k, k + 1, tol);
		// std::cout << "#####################" << "\n";
		// printSubMatrix(B, p, q, p, q);
		// std::cout << "#####################" << "\n";

		// std::cout << "U\n";
		// std::cout << c << " " << s << " " << k << " " << k+1 << "\n";
		// printMatrix(uT);
		// std::cout << "IS\n";
		// printMatrix(uT*B*v);
		y = B[k][k + 1]; z = B[k][k + 2];
	
		givens(y, z, &c, &s, tol);
		rightGivens(B, c, s, k + 1, k + 2, tol);
		leftGivens(v, c, s, k + 1, k + 2, tol);
	
		// std::cout << "#####################" << "\n";
		// printSubMatrix(B, p, q, p, q);
		// std::cout << "#####################" << "\n";

	}

	// std::cout << "#####################" << "\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "##########AAA###########" << "\n";

	y = B[q - 2][q - 2]; z = B[q - 1][q - 2];
	givens(y, z, &c, &s, tol);
	// std::cout << y << " " << z << "\n";
	// std::cout << c << " " << s << "\n";
	leftGivens(B, c, s, q - 2, q - 1, tol);
	rightGivens(uT, c, s, q - 2, q - 1, tol);
	// printSubMatrix(B, p, q, p, q);

	// std::cout << "##########AAA###########" << "\n";

	// std::cout << "B Matrix:\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "\n";

}

void clean(Matrix& B, f32 tol)
{
	i32 i, n;
	n = B.columns();
	for(i=0; i < n; i++)
	{
		if(fabs(B[i][i]) < tol)
		{
			B[i][i] = 0.0;
		}
	}
	for(i=0; i < n - 1; i++)
	{
		if(fabs(B[i][i+1]) < tol * (fabs(B[i][i]) + fabs(B[i + 1][i + 1])))
		{
			B[i][i + 1] = 0.0;
		}
	}
}

void findLimits(Matrix& B, i32* p, i32 *q, f32 tol)
{
	i32 i, n;
	
	n = B.columns();

	*q = -1;
	*p = -1;

	// std::cout << "************ findLimits ************\n";
	
	for(i = n - 1; i >=1; i--)
	{
		// Pick largest q such that B33 is diagonal
		// std::cout << "(" << i - 1 << ", " << B[i][i] << "), ";
		
		if(*q == -1 && B[i - 1][i])
		{
			*p = *q = i - 1;
		}
		else if(*q != -1)
		{
			if(B[i - 1][i])
			{
				*p = i - 1;
			}
			else
			{
				break;
			}
		}
	}

	// std::cout << "\n";

	if(*q == -1)
	{
		*p = 0;
		*q = 0;
	}
	else
	{
		*p = *p;
		*q = *q + 1;
	}

	// std::cout << *p << " " << *q << "\n";
}

void clearDiag(f32* b_diag, f32* b_sdiag, i32 m, i32 n, f32 tol)
{
	i32 i;
	for(i = 0; i < n - 1; i++)
	{
		if(fabs(b_sdiag[i+1]) <= tol*(fabs(b_diag[i]) + fabs(b_diag[i+1])))
		{
			b_sdiag[i+1] = 0.0;
		}
	}
	
	for(i = 0; i < n; i++)
	{
		if(fabs(b_diag[i]) <= tol)
			b_diag[i] = 0.0;
		if(fabs(b_sdiag[i]) <= tol)
			b_sdiag[i] = 0.0;
	}
}

f32 getDiagMat(f32** diags, i32 l, i32 c)
{
	if((l + 2 - c >= 0) && (l + 2 - c < 5))
	{
		return diags[l + 2 - c][c];
	}

	return 0;
}


void walkBlemishOutRight(Matrix& B, i32 p, i32 q, i32 row, i32 start_col, Matrix& uT, Matrix& v, f32 tol)
{
	Matrix B_ = Matrix(B);

	i32 col, n;
	f32 c, s, brr, brc, old;

	n = B.columns();

	for(col = start_col; col < q; col++)
	{
		// b[row][row]
		// std::cout << "out right:\n";
		// std::cout << "[" << col << ", " << col << "] ";
		// std::cout << brr << " *** ";
		// std::cout << "[" << row << ", " << col << "] ";
		// std::cout << brc << "\n";

		// std::cout << "B[col][col] << " " << B[row][col]" << "\n";
		givens(B[col][col], B[row][col], &c, &s, tol);
		// std::cout << std::scientific << B[col][col] << " " << B[row][col] << "\n";
		old = B[row][col];

		B[row][col] = 0;                            
		B[col][col] = B[col][col] * c - old * s;
		// if(fabs(B[col][col]) <= tol) B[col][col] = 0.0;

		if(col < n - 1)
		{
			B[row][col + 1] = s * B[col][col + 1];
			// if(fabs(B[row][col + 1]) <= tol) B[row][col + 1] = 0.0;
			B[col][col + 1] = c * B[col][col + 1];
			// if(fabs(B[col][col + 1]) <= tol) B[col][col + 1] = 0.0;
		}

		rightGivens(uT, c, s, col, row, tol);
		
		// std::cout << "U in blemish\n";
		// std::cout << std::scientific << c << " " << s << " " << col << " " << row << "\n";
		// printMatrix(uT);
		// std::cout << "#####################*****************\n";
	}

	// std::cout << "IS:\n";
	// printMatrix(uT*B*v);

	// std::cout << "\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "\n";
}

// start row is kk-1 and col kk sup diag od k
void walkBlemishOutUp(Matrix& B, i32 p, i32 q, i32 start_row, i32 col, Matrix& uT, Matrix& v, f32 tol)
{
	i32 row;
	f32 c, s, old;
	// std::cout << "walkBlemishOutUp:\n";
	// std::cout << start_row << " " << p << "\n";

	for(row = start_row; row >= p; row--)
	{
		// std::cout << "left up:\n";
		// std::cout << "[" << row << ", " << row << "] ";
		// std::cout << B[row][row] << " *** ";
		// std::cout << "[" << row << ", " << col << "] ";
		// std::cout << B[row][col] << "\n";

		givens(B[row][row], B[row][col], &c, &s, tol);

		old = B[row][col];
	
		B[row][col] = 0;    
	                        
		B[row][row] = B[row][row] * c - old * s;

		// if(fabs(B[row][row]) <= tol) B[row][row] = 0.0;

		if(row > 0)
		{
			B[row-1][col] = s * B[row - 1][row];
			// if(fabs(B[row-1][col]) <= tol) B[row-1][col] = 0.0;
		
			B[row-1][row] = c * B[row-1][row];
			// if(fabs(B[row-1][row]) <= tol) B[row-1][row] = 0.0;
		}

		leftGivens(v, c, s, row, col, tol);

		// std::cout << "V in blemish\n";
		// std::cout << c << " " << s << " " << row << " " << col << "\n";
	}

	// std::cout << "IS:\n";
	// printMatrix(uT*B*v);
	// std::cout << "\n";
	// printSubMatrix(B, p, q, p, q);
	// std::cout << "\n";
	// std::cout << "***************************\n";

}


// GvL pg 454
bool doZeroDiag(Matrix& B, i32 p, i32 q, Matrix& uT, Matrix& v, f32 tol)
{
	i32 i = 0;

	bool zeroed = false;

	// std::cout << "doZeroDiag" "\n";
	// std::cout << p << " " << q << "\n";

	// printSubMatrix(B, p, q, p, q);
	// Get zeros idx in diag
	for(i = p; i < q; i++)
	{

		if(fabs(B[i][i]) <= tol)
		{
			// std::cout << "zeroidx: " << i << "\n";
			zeroed = true;
			
			if(i < q - 1)
			{
				walkBlemishOutRight(B, p, q, i, i+1, uT, v, tol);
			}
			else if(i == q - 1)
			{
				walkBlemishOutUp(B, p, q,  i-1, i, uT, v, tol);
			}
		}
	}

	return zeroed;
}


// void _golubKahanSVD(f32* b_diag, f32* b_sdiag, i32 m, i32 n, Matrix &uT, Matrix &v, f32 tol)
// {
// 	i32 i, j, q = 0, p = 0;
// 	uT = transpose(uT);

// 	while(q != n)
// 	{
// 		// Remove elements that are less that the tolerance
// 		clearDiag(b_diag, b_sdiag, m, n, tol);
		
// 		// Find the largest q and the smallest p such that if
// 		// B = [B11, 0, 0; 0, B22, 0; 0, 0, B33] then B33 is diagonal
// 		// and B22 has nonzero superdiagonal
// 		findLimits(b_diag, b_sdiag, m, n, &p, &q, tol);

// 		// B11 = B[0:p-1, 0:p-1]
// 		// B22 = B[p:n - q - 1, p:n - q - 1]
// 		// B33 = B[n - q:n, n - q:n]

// 		if(q < n)
// 		{
// 			// if any entry in B22 is zero, then zero
// 			// the superdiagonal entry in the same row.
// 			bool zeroed = false;

// 			for(i=0; i < n - p - q; i++)
// 			{
// 				if(i < n - 1 && fabs(b_diag[i]) <= tol)
// 				{
// 					b_sdiag[i+1] = 0.0;
// 					zeroed = true;
// 				}
// 			}

// 			if(!zeroed)
// 			{
// 				f32* b22_diag  = new f32[n - p - q];
// 				f32* b22_sdiag = new f32[n - p - q];

// 				Matrix uT_ = identity(n - p - q, n - p - q);
// 				Matrix v_  = identity(n - p - q, n - p - q);

// 				for(i = 0; i < n - p - q; i++)
// 				{
// 					b22_diag[i]  = b_diag[i + p];
// 					b22_sdiag[i] = b_sdiag[i + p];
// 				}

// 				Matrix Iu = identity(p + (n - p - q) + (q + m - n), p + (n - p - q) + (q + m - n));
// 				Matrix Iv = identity(p + (n - p - q) + q, p + (n - p - q) + q);
// 				for(i=p; i < n - p - q; i++)
// 				{
// 					for(j=p; j < n - p - q; j++)
// 					{
// 						Iu[i][j] = uT_[i - p][j - p];
// 					}
// 				}

// 				for(i=p; i < n - p - q; i++)
// 				{
// 					for(j=p; j < n - p - q; j++)
// 					{
// 						Iv[i][j] = v_[i][j];
// 					}
// 				}

// 				uT = Iu*uT; v  = v*Iv;


// 				for(i = 0; i < n - p - q; i++)
// 				{
// 					b_sdiag[i + p] = b22_sdiag[i];
// 					b_diag[i + p]  = b22_diag[i];
// 				}

// 				delete[] b22_diag;
// 				delete[] b22_sdiag;
// 			}
// 		}
// 	}

// 	uT  = transpose(uT);

// 	for(i=0; i<n; i++)
// 	{
// 		if(b_diag[i] < 0)
// 		{
// 			b_diag[i] = -b_diag[i];

// 			if(i%2 == 1)
// 			{
// 				for(j=0; j<n; j++)
// 				{
// 					v[j][i] = -v[j][i];
// 				}
// 			}
// 			else
// 			{
// 				for(j=0; j<m; j++)
// 				{
// 					uT[j][i] = -uT[j][i];
// 				}
// 			} 
// 		}
// 	}

// 	if(uT[0][0] > 0)
// 	{
// 		for(j=0; j<m; j++)
// 		{
// 			uT[j][0] = -uT[j][0];
// 		}
// 		for(j=0; j<m; j++)
// 		{
// 			v[0][j] = -v[0][j];
// 		}
// 	}

// 	sortUV(b_diag, n, uT, v);

// }


void golubKahanSVD(f32* b_diag, f32* b_sdiag, i32 m, i32 n, Matrix &uT, Matrix &v, f32 tol)
{
	i32 i, j, q = 0, p = 0;
	
	Matrix B = bidiag(b_diag, b_sdiag, m, n);
	
	// std::cout << "U\n";
	// printMatrix(uT);
	// std::cout << "\n";
	// std::cout << "V\n";
	// printMatrix(v);
	// std::cout << "CLEAN\n";
	v = transpose(v);

	// std::cout << "**************B*************\n";
	// printMatrix(B);
	// std::cout << "**************B*************\n";
	clean(B, /* tol */ tol);
	// printMatrix(B);
	// std::cout << "CLEAN\n";
	// std::cout << "**************B*************\n";
	// printMatrix(B);
	// std::cout << "**************B*************\n";
	findLimits(B, &p, &q, tol);
	// printSubMatrix(B, p, q+1, p, q+1);
	// std::cout << "************** START *************\n";

	while(q - p)
	{

		// std::cout << p << " " << q << "\n";
		// std::cout << "**************B*************\n";
		// printMatrix(B);
		// std::cout << "**************B*************\n";
		bool zeroed = doZeroDiag(B, p, q + 1, uT, v, tol);
		if(!zeroed)
		{
			// std::cout << "Golub step\n";
			golubKahanStep(B, p, q + 1, uT, v, tol);
		}
		else
		{
			// std::cout << "zeroes\n";
		}

		clean(B, /* tol */ tol);
		findLimits(B, &p, &q, tol);

		// std::cout << "************* END ITERATION **************\n" << "\n";
	}

	// v = transpose(v);

	for(i=0; i<n; i++)
	{
		b_diag[i] = B[i][i];
	}

	sortUV(b_diag, n, uT, v);

	// std::cout << "Bdiag\n";

	// printMatrix(diag(b_diag, m, n));

	// std::cout << "U\n";

	// printMatrix(uT);
	// std::cout << "V\n";
	// printMatrix(v);
	
	// std::cout << "U*S*VT\n";
	// printMatrix(uT*diag(b_diag, m, n)*v);
	// std::cout << "\n";
}


// f32 sigmaMin(f32* s, f32* e, i32 m, i32 n)
// {
// 	i32 j;
// 	f32 *l, *u, b_inf, b_one, sigma;

// 	l = (f32*)malloc(sizeof(f32)*n);
// 	u = (f32*)malloc(sizeof(f32)*n);

// 	l[n-1] = fabs(s[n-1]);

// 	for(j = n-2; j >= 0; j--)
// 		l[j] = fabs(s[j]) * (l[j+1] / (l[j+1] + fabs(e[j])));

// 	u[0] = fabs(s[0]);

// 	for(j=0; j<n-1; j++)
// 		u[j+1] = fabs(s[j+1]) * (u[j] / (u[j] + fabs(e[j])));

// 	b_inf = l[0];
// 	b_one = u[0];

// 	for(j=1; j<n; j++)
// 		b_inf = std::min(b_inf, l[j]);

// 	for(j=1; j<n; j++)
// 		b_one = std::min(b_one, u[j]);

// 	sigma = std::min(b_one, b_inf);

// 	delete[] l;
// 	delete[] u;
// 	std::cout << "sigma: " << sigma << "\n";
// 	return sigma;
// }



}
